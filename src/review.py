"""Read-only inspection and opt-in review edits for saved namiyomi projects.

Used by the MCP server, but has no MCP dependency so it can be tested directly.
No assembly or consensus logic lives here. Review decisions go into the same per-contig store the
in-app AI review writes (``Contig.ai_review`` / ``Contig.ai_notes``), so the app shows them like any
other review result: a decided column, or a red "unresolved" flag. End trims use ``Contig.overrides``
(like deleting a consensus base by hand). Everything goes to a NEW project file, never over the input.
"""
import io
import shutil
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from types import SimpleNamespace

from batch import CLEAN, NEEDS_REVIEW, SINGLE_READ, triage
from core import COMPLEMENT, Project, atomic_text, base_text, export_markers
from sequence_find import consensus_coordinates, find_matches
from trace_plot import trace_samples

# Alignment columns are numbered from 1 everywhere in the tools, reports and tool arguments, the
# same as the app's column ruler and review.tsv. Internally they stay 0-based.
def _one(column):
    return None if column is None else column + 1


def _zero(value, what='column'):
    try:
        column = int(value)
    except (TypeError, ValueError):
        raise ValueError(f'{what} must be a whole number (columns are numbered from 1).') from None
    return column - 1


BASE_COLORS = {'A': '#008234', 'C': '#1559c5', 'G': '#161616', 'T': '#d1272c'}


def load_project(path):
    return Project.load(Path(path).expanduser())


def resolve_contig(project, ref):
    """Find a contig by display name (as in summary.json) or by its internal key."""
    for key, contig in project.contigs.items():
        if ref in (contig.sequence_name, key):
            return key, contig
    names = ', '.join(c.sequence_name for c in project.contigs.values())
    raise ValueError(f'No contig named {ref!r}. Contigs: {names}')


def list_contigs(project, category=None):
    rows = triage(project)
    if category:
        if category not in (CLEAN, NEEDS_REVIEW, SINGLE_READ):
            raise ValueError(f'Unknown filter {category!r}; use clean, needs_review or single_read.')
        rows = [r for r in rows if r['category'] == category]
    return rows


def contig_details(project, ref):
    """Flagged columns (alignment columns, numbered from 1 like the app) with every source read's call and Phred."""
    key, contig = resolve_contig(project, ref)
    row = next(r for r in triage(project) if r['key'] == key)
    sources = []
    for index, source in enumerate(contig.sources):
        read = project.reads[source.read_id]
        sources.append(dict(index=index, file=Path(read.path).name, reverse_complemented=source.reverse,
                            trim_start=read.start + 1, trim_end=read.end))
    flagged = []
    for col in contig.issues():
        reads = []
        for index, source in enumerate(contig.sources):
            b = source.mapping[col]
            reads.append(dict(source=index, base=None if b is None else b.base,
                              phred=None if b is None else b.quality,
                              raw_position=None if b is None or b.raw is None else b.raw + 1,
                              peak_heights=_peak_heights(project.traces[source.read_id], b, source.reverse)))
        column = contig.columns[col]
        flagged.append(dict(column=_one(col), consensus=contig.call(col), automatic_call=column.call,
                            reason=column.reason, manual_override=str(col) in contig.overrides, reads=reads))
    return dict(name=contig.sequence_name, key=key, category=row['category'], status=contig.status,
                overlap=contig.overlap, identity=row['identity'], columns=len(contig.columns),
                min_quality=project.settings['min_quality'], sources=sources, flagged=flagged)


def column_qualities(project, ref, start=None, end=None):
    """Every alignment column in [start, end] (numbered from 1, inclusive): consensus call and each read's base and Phred.

    Use it to choose a trim point, e.g. where a rolling mean of the Phred scores falls below a threshold.
    """
    _, contig = resolve_contig(project, ref)
    total = len(contig.columns)
    start = 0 if start is None else _zero(start, 'start')
    end = total - 1 if end is None else _zero(end, 'end')
    if not 0 <= start <= end < total:
        raise ValueError(f'Need 1 <= start <= end <= {total}.')
    rows = []
    for col in range(start, end + 1):
        reads = [None if b is None else [b.base, b.quality] for b in (s.mapping[col] for s in contig.sources)]
        rows.append(dict(column=_one(col), consensus=contig.call(col), reads=reads))
    return dict(name=contig.sequence_name, columns=total, min_quality=project.settings['min_quality'],
                sources=[Path(project.reads[s.read_id].path).name for s in contig.sources],
                note='reads: [base, phred] per source in that order, null if the read does not cover the column',
                rows=rows)


def _rolling_keep(scores, window, threshold, limit):
    """First index (scanning from 0, never past ``limit``) whose window mean reaches ``threshold``.

    ``scores`` holds None for columns left out of the mean (isolated Ns, bases without a Phred).
    """
    for i in range(limit):
        vals = [q for q in scores[i:i + window] if q is not None]
        if not vals or sum(vals) / len(vals) >= threshold:
            return i
    return limit


def _n_run_cut(calls, run, cols):
    """Walk ``cols`` (inner to outer) and return the innermost column of the first run of ``run`` Ns.

    Gap columns ('' or '-') are skipped, so they neither break nor extend a run. None if there is no such run.
    """
    count, first = 0, None
    for col in cols:
        if calls[col] in ('', '-'):
            continue
        if calls[col] == 'N':
            count, first = count + 1, col if count == 0 else first
            if count >= run:
                return first
        else:
            count = 0
    return None


class _Ends:
    """Per-column data for one contig, shared by the end-trimming helpers."""

    def __init__(self, contig, keep_isolated_n=True):
        n = len(contig.columns)
        self.n, self.contig = n, contig
        self.calls = [contig.call(i) for i in range(n)]
        self.scores, self.phred, self.covered = [], [], []
        for col in range(n):
            bases = [s.mapping[col] for s in contig.sources if s.mapping[col] is not None]
            quals = [b.quality for b in bases if b.quality is not None]
            self.covered.append(len(bases))
            self.phred.append(max(quals) if quals else None)
            isolated_n = (keep_isolated_n and self.calls[col] == 'N'
                          and (col == 0 or self.calls[col - 1] != 'N')
                          and (col == n - 1 or self.calls[col + 1] != 'N'))
            self.scores.append(None if isolated_n else self.phred[-1])
        multi = [i for i, c in enumerate(self.covered) if c >= 2]
        self.overlap = (multi[0], multi[-1]) if multi else None

    def bases(self, cols):
        return [self.calls[i] for i in cols if self.calls[i] not in ('', '-')]

    def bounds(self, window=10, min_mean_phred=20, single_read_ends_only=True, max_n_run=None):
        """(keep_from, keep_to) under the rules; raises ValueError if nothing would be kept."""
        _check_rule(window, min_mean_phred, max_n_run)
        n, restrict = self.n, single_read_ends_only and self.overlap
        left_limit = self.overlap[0] if restrict else n
        right_limit = n - 1 - self.overlap[1] if restrict else n
        keep_from, keep_to = 0, n - 1
        if min_mean_phred is not None:
            keep_from = _rolling_keep(self.scores, window, float(min_mean_phred), left_limit)
            keep_to = n - 1 - _rolling_keep(self.scores[::-1], window, float(min_mean_phred), right_limit)
        if max_n_run:
            # Walk outward from the reliable middle: the overlap if restricted, else the contig's centre.
            inner_left, inner_right = (self.overlap[0] - 1, self.overlap[1] + 1) if restrict else (n // 2 - 1, n // 2)
            cut = _n_run_cut(self.calls, int(max_n_run), range(inner_right, n))
            if cut is not None:
                keep_to = min(keep_to, cut - 1)
            cut = _n_run_cut(self.calls, int(max_n_run), range(inner_left, -1, -1))
            if cut is not None:
                keep_from = max(keep_from, cut + 1)
        if keep_from > keep_to:
            raise ValueError('These rules would remove the whole contig.')
        return keep_from, keep_to

    def outcome(self, keep_from, keep_to):
        kept = self.bases(range(keep_from, keep_to + 1))
        return dict(keep_from=_one(keep_from), keep_to=_one(keep_to), length=len(kept), n_count=kept.count('N'),
                    removed_left=len(self.bases(range(keep_from))),
                    removed_right=len(self.bases(range(keep_to + 1, self.n))))

    def tail(self, cols):
        cols = list(cols)
        bases = self.bases(cols)
        quals = [self.phred[i] for i in cols if self.phred[i] is not None and self.calls[i] not in ('', '-')]
        n_cols = [i for i in cols if self.calls[i] == 'N']
        longest = run = 0
        for i in cols:
            if self.calls[i] in ('', '-'):
                continue
            run = run + 1 if self.calls[i] == 'N' else 0
            longest = max(longest, run)
        return dict(columns=[_one(cols[0]), _one(cols[-1])] if cols else None, bases=len(bases),
                    mean_phred=round(sum(quals) / len(quals), 1) if quals else None,
                    n_count=len(n_cols), longest_n_run=longest, n_columns=[_one(i) for i in n_cols])


def _check_rule(window, min_mean_phred, max_n_run):
    if int(window) < 1:
        raise ValueError('window must be at least 1.')
    if max_n_run is not None and int(max_n_run) < 1:
        raise ValueError('max_n_run must be at least 1.')
    if min_mean_phred is None and not max_n_run:
        raise ValueError('Give min_mean_phred, max_n_run or both.')


def _rule_text(window, min_mean_phred, single_read_ends_only, keep_isolated_n, max_n_run):
    parts = []
    if min_mean_phred is not None:
        parts.append(f'rolling mean Phred over {window} columns below {min_mean_phred:g}'
                     + (', isolated Ns left out of the mean' if keep_isolated_n else ''))
    if max_n_run:
        parts.append(f'cut before the first run of {max_n_run} or more Ns')
    scope = 'single-read ends only' if single_read_ends_only else 'both ends'
    return '; '.join(parts) + f' ({scope}).' if parts else 'no rule'


def _contigs(project, ref):
    return [resolve_contig(project, ref)] if ref else list(project.contigs.items())


def suggest_quality_trim(project, ref=None, window=10, min_mean_phred=20, single_read_ends_only=True,
                         keep_isolated_n=True, max_n_run=None):
    """Dry run of end trimming of contigs. Nothing is changed.

    Two rules, either or both (the stricter cut wins at each end):

    - Rolling mean (``min_mean_phred``; None turns it off): column quality is the highest Phred among the
      reads covering it; from each end inward, columns are trimmed while the mean over the next
      ``window`` columns is below the threshold. With ``keep_isolated_n``, an N flanked by non-N calls is
      left out of the mean, so it neither causes nor blocks a cut.
    - N runs (``max_n_run``): walking outward from the region two reads cover, each end is cut just before
      the first run of ``max_n_run`` or more consecutive Ns.

    With ``single_read_ends_only``, cuts never reach into the region two reads cover (contigs with one read
    are trimmed from their centre outward). Each result carries a ready ``trim`` edit for
    :func:`apply_review_edits` (``None`` when nothing would be trimmed).
    """
    window = int(window)
    _check_rule(window, min_mean_phred, max_n_run)
    rule = _rule_text(window, min_mean_phred, single_read_ends_only, keep_isolated_n, max_n_run)
    results = []
    for _, contig in _contigs(project, ref):
        if contig.stale:
            results.append(dict(name=contig.sequence_name, error='Inputs have changed; rebuild the contig first.'))
            continue
        ends = _Ends(contig, keep_isolated_n)
        try:
            keep_from, keep_to = ends.bounds(window, min_mean_phred, single_read_ends_only, max_n_run)
        except ValueError as exc:
            results.append(dict(name=contig.sequence_name, columns=ends.n, error=str(exc)))
            continue
        edit = None
        if (keep_from, keep_to) != (0, ends.n - 1):
            edit = dict(contig=contig.sequence_name, action='trim', keep_from=_one(keep_from), keep_to=_one(keep_to),
                        reason=rule[0].upper() + rule[1:])
        result = dict(name=contig.sequence_name, columns=ends.n,
                      overlap=[_one(c) for c in ends.overlap] if ends.overlap else None, **ends.outcome(keep_from, keep_to))
        result['isolated_n_kept'] = sum(1 for i in range(keep_from, keep_to + 1)
                                        if ends.calls[i] == 'N' and ends.scores[i] is None) if keep_isolated_n else 0
        result['edit'] = edit
        results.append(result)
    return results


DEFAULT_TRIM_OPTIONS = (
    dict(label='keep everything'),
    dict(label='rolling mean Q20 over 10', window=10, min_mean_phred=20),
    dict(label='rolling mean Q30 over 20', window=20, min_mean_phred=30),
    dict(label='cut before a run of >=4 N', min_mean_phred=None, max_n_run=4),
    dict(label='cut before a run of >=3 N', min_mean_phred=None, max_n_run=3),
    dict(label='cut before any N', min_mean_phred=None, max_n_run=1),
    dict(label='region covered by two reads', overlap_only=True),
)


def describe_contig_ends(project, ref=None, options=None):
    """What the ends of each contig look like and what candidate trimming rules would do. Nothing is changed.

    For each contig: how the reads are laid out, the region two reads cover, each single-read end
    (length, mean Phred, N count, longest N run, N columns), and for each option the kept range, length
    and N count. ``options`` is a list of dicts with ``label`` plus any of the
    :func:`suggest_quality_trim` arguments (``window``, ``min_mean_phred``, ``max_n_run``,
    ``single_read_ends_only``, ``keep_isolated_n``), or ``overlap_only: true``; an option with no rule
    keeps everything. The default compares common rules.
    """
    options = list(options or DEFAULT_TRIM_OPTIONS)
    rows = {r['key']: r for r in triage(project)}
    results = []
    for key, contig in _contigs(project, ref):
        ends = _Ends(contig)
        reads = []
        for source in contig.sources:
            read, trace = project.reads[source.read_id], project.traces[source.read_id]
            cols = [i for i, b in enumerate(source.mapping) if b is not None]
            kept = [q for q in trace.qualities[read.start:read.end] if q is not None]
            reads.append(dict(file=Path(read.path).name, reverse_complemented=source.reverse,
                              columns=[_one(cols[0]), _one(cols[-1])] if cols else None, raw_length=len(trace.sequence),
                              trimmed_length=read.end - read.start,
                              median_phred=sorted(kept)[len(kept) // 2] if kept else None))
        if ends.overlap:
            left, right = range(0, ends.overlap[0]), range(ends.overlap[1] + 1, ends.n)
        else:
            left = right = range(0)
        compared = []
        for option in options:
            label = option.get('label') or 'option'
            try:
                if option.get('overlap_only'):
                    if not ends.overlap:
                        continue
                    bounds = ends.overlap
                elif not any(k in option for k in ('min_mean_phred', 'max_n_run')):
                    bounds = (0, ends.n - 1)
                else:
                    kw = {k: option[k] for k in ('window', 'min_mean_phred', 'single_read_ends_only', 'max_n_run')
                          if k in option}
                    bounds = _Ends(contig, option.get('keep_isolated_n', True)).bounds(**kw)
                compared.append(dict(label=label, **ends.outcome(*bounds)))
            except (TypeError, ValueError) as exc:
                compared.append(dict(label=label, error=str(exc)))
        results.append(dict(name=contig.sequence_name, category=rows[key]['category'], columns=ends.n,
                            overlap=[_one(c) for c in ends.overlap] if ends.overlap else None, reads=reads,
                            left_end=ends.tail(left) if ends.overlap else None,
                            right_end=ends.tail(right) if ends.overlap else None,
                            n_total=ends.calls.count('N'), options=compared))
    return results


def _peak_heights(trace, base, reverse):
    """Raw height of each channel at the base's peak, named in contig orientation (None if no peak)."""
    if base is None or base.raw is None:
        return None
    x = trace.peaks[base.raw]
    return {COMPLEMENT[ch] if reverse else ch: values[x] for ch, values in sorted(trace.channels.items())}


def find_sequence(project, ref, query, identity=100):
    """Matches in the exported consensus: start..end are positions from 1, inclusive; *_column are alignment columns."""
    _, contig = resolve_contig(project, ref)
    sequence, columns = consensus_coordinates(contig)
    return [dict(start=m.start + 1, end=m.end, strand=m.strand, identity=m.identity,
                 first_column=_one(columns[m.start]), last_column=_one(columns[m.end - 1]))
            for m in find_matches(sequence, query, identity)]


def render_trace(project, ref, column, window=15):
    """PNG (bytes) of every source read's chromatogram around one alignment column (numbered from 1)."""
    from matplotlib.figure import Figure
    key, contig = resolve_contig(project, ref)
    total = len(contig.columns)
    column = _zero(column)
    if not 0 <= column < total:
        raise ValueError(f'Column must be between 1 and {total}.')
    window = max(3, min(int(window), 60))
    left, right = max(0, column - window), min(total, column + window + 1)
    fig = Figure(figsize=(10, 2.6 * len(contig.sources)), dpi=100)
    axes = fig.subplots(len(contig.sources), 1, sharex=True, squeeze=False)[:, 0]
    minq = project.settings['min_quality']
    for ax, source in zip(axes, contig.sources):
        read, trace = project.reads[source.read_id], project.traces[source.read_id]
        samples = trace_samples(trace, source.mapping, left, right, 1.0, source.reverse)
        xs = [(x - 82) + left for x, _ in samples]
        peak = max((trace.channels[c][s] for c in trace.channels for _, s in samples), default=1) or 1
        for channel, values in trace.channels.items():
            shown = COMPLEMENT[channel] if source.reverse else channel
            ax.plot(xs, [values[s] / peak for _, s in samples], color=BASE_COLORS[shown], lw=1.1, label=shown)
        ax.axvspan(column - .5, column + .5, color='#f0c040', alpha=.3)
        ax.axhline(0, color='#ddd', lw=.8)
        for col in range(left, right):
            b = source.mapping[col]
            if b is None:
                ax.text(col, 1.12, '-', ha='center', color='#777', fontsize=9)
                continue
            ax.text(col, 1.12, b.base or '-', ha='center', fontsize=10, fontweight='bold',
                    color=BASE_COLORS.get(b.base, '#555'))
            if b.quality is not None:
                ax.text(col, 1.0, str(b.quality), ha='center', fontsize=6,
                        color='#222' if b.quality >= minq else '#b00020')
        ax.set_ylim(-.05, 1.25)
        ax.set_yticks([])
        ax.set_xlim(left - .5, right - .5)
        strand = ' (reverse-complemented)' if source.reverse else ''
        ax.set_title(f'{Path(read.path).name}{strand}', loc='left', fontsize=9)
    axes[-1].set_xlabel(f"{contig.sequence_name}: alignment column {column + 1} (highlighted). "
                        f"Top row = base call, then Phred (red < Q{minq}). Consensus = {contig.call(column)!r}")
    fig.tight_layout()
    buf = io.BytesIO()
    fig.savefig(buf, format='png')
    return buf.getvalue()


NOTE_LIMIT = 200   # Project.load rejects longer review notes


def _note(reviewer, reason):
    text = ' '.join(f'{reviewer} (MCP): {reason}'.split())
    text = ''.join(c for c in text if c >= ' ')
    return text if len(text) <= NOTE_LIMIT else text[:NOTE_LIMIT - 1] + '…'


def apply_review_edits(project_path, output_path, edits, report_path=None, reviewer='Claude', unresolved=None):
    """Apply review verdicts to a COPY of the project and log each with its reason.

    Each edit: ``{contig, column, action: 'set'|'confirm', value (for set), reason}``, or
    ``{contig, action: 'trim', keep_from, keep_to, reason}``.
    ``set`` records a review decision for the column (A, C, G or T; ``-`` removes the column);
    ``confirm`` marks the flag reviewed without changing the call; ``trim`` deletes the consensus columns outside
    ``keep_from..keep_to`` (inclusive, like deleting bases in the GUI's consensus row).
    All edits are validated before any is applied.
    ``unresolved`` optionally gives ``{contig, column, reason}`` for flags deliberately left alone;
    every flag still open after the edits is listed in the report, with its reason if given.
    """
    source, output = Path(project_path).expanduser().resolve(), Path(output_path).expanduser().resolve()
    if source == output:
        raise ValueError('The reviewed project must be saved to a new file; the input is never overwritten.')
    edits, unresolved = list(edits or []), list(unresolved or [])
    if not edits and not unresolved:
        raise ValueError('No edits given.')
    project = load_project(source)
    plan = []
    for n, edit in enumerate(edits, start=1):
        try:
            key, contig = resolve_contig(project, edit['contig'])
            action = edit.get('action', 'set')
            reason = str(edit.get('reason', '')).strip()
            last = len(contig.columns) - 1
            total = last + 1
            if contig.stale:
                raise ValueError('Inputs have changed; rebuild the contig first.')
            if action not in ('set', 'confirm', 'trim'):
                raise ValueError("action must be 'set', 'confirm' or 'trim'.")
            if not reason:
                raise ValueError('A reason is required so the change can be audited.')
            value = None
            if action != 'trim' and str(_zero(edit['column'])) in contig.overrides:
                raise ValueError('That column was edited by hand; leave it alone.')
            if action == 'trim':
                col = (_zero(edit['keep_from'], 'keep_from'), _zero(edit['keep_to'], 'keep_to'))
                if not 0 <= col[0] <= col[1] <= last:
                    raise ValueError(f'Need 1 <= keep_from <= keep_to <= {total}.')
                if col == (0, last):
                    raise ValueError('keep_from..keep_to covers the whole contig; nothing to trim.')
            else:
                col = _zero(edit['column'])
                if not 0 <= col <= last:
                    raise ValueError(f'Column {col + 1} is outside 1..{total}.')
            if action == 'set':
                value = str(edit.get('value', '')).upper().strip()
                if value not in ('A', 'C', 'G', 'T', '-'):
                    raise ValueError("value must be A, C, G, T, or '-' to remove the column.")
            plan.append((key, contig, col, action, value, reason))
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError(f'Edit {n}: {exc.args[0] if isinstance(exc, KeyError) else exc}') from exc
    notes, unresolved_marks = {}, []
    for n, item in enumerate(unresolved, start=1):
        try:
            _, contig = resolve_contig(project, item['contig'])
            col = _zero(item['column'])
            if not 0 <= col < len(contig.columns):
                raise ValueError(f'Column {col + 1} is outside 1..{len(contig.columns)}.')
            reason = str(item.get('reason', '')).strip()
            if str(col) not in contig.overrides:
                unresolved_marks.append((contig, col, reason))
            notes[(contig.sequence_name, col)] = reason
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError(f'Unresolved item {n}: {exc.args[0] if isinstance(exc, KeyError) else exc}') from exc
    lines = []
    for key, contig, col, action, value, reason in plan:
        if action == 'trim':
            keep_from, keep_to = col
            removed = [i for i in range(len(contig.columns)) if not keep_from <= i <= keep_to]
            bases = sum(contig.call(i) not in ('', '-') for i in removed)
            for i in removed:
                contig.overrides[str(i)] = ''          # same as deleting a consensus base in the GUI
                contig.black_marks.pop(str(i), None)
            cut = ([f'1-{keep_from}'] if keep_from else []) + \
                  ([f'{keep_to + 2}-{len(contig.columns)}'] if keep_to < len(contig.columns) - 1 else [])
            lines.append(f'- **{contig.sequence_name}**: trimmed to columns {keep_from + 1}-{keep_to + 1} '
                         f'({bases} bases removed from columns {" and ".join(cut)}). {reason}')
            continue
        before = contig.call(col)
        if action == 'set':
            contig.ai_review[str(col)] = '' if value == '-' else value
            contig.ai_notes[str(col)] = _note(reviewer, reason)
            contig.black_marks.pop(str(col), None)
            shown = 'no base (column removed)' if value == '-' else value
            lines.append(f'- **{contig.sequence_name}** column {col + 1}: {before} -> **{shown}**. {reason}')
        else:
            contig.black_marks[str(col)] = False    # same as the GUI's "Mark reviewed"
            contig.ai_notes.pop(str(col), None)
            lines.append(f'- **{contig.sequence_name}** column {col + 1}: confirmed {before}. {reason}')
    for contig, col, reason in unresolved_marks:       # shown as a red flag in the app
        if contig.ai_review.get(str(col), 'N') == 'N':
            contig.ai_review[str(col)] = 'N'
            contig.ai_notes[str(col)] = _note(reviewer, reason or 'left for the user')
    project.dirty = True
    project.save(output)
    remaining = [r for r in triage(project) if r['category'] == NEEDS_REVIEW]
    open_lines = []
    for r in remaining:
        for col in r['issue_columns']:
            note = notes.get((r['name'], col))
            open_lines.append(f'- **{r["name"]}** column {col + 1}' + (f': {note}' if note else ''))
    body = '\n'.join(lines) if lines else '_No edits._'
    body += '\n\n### Left unresolved (check these in the app)\n\n' + ('\n'.join(open_lines) or '_None._')
    report = Path(report_path).expanduser() if report_path else output.with_name('review_report.md')
    stamp = datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC')
    old = report.read_text(encoding='utf-8') if report.exists() else '# Review report\n'
    atomic_text(report, old.rstrip('\n') + f'\n\n## {reviewer}, {stamp}\n'
                f'Source project: `{source.name}` (unchanged). Reviewed copy: `{output.name}`.\n\n'
                + body + '\n')
    return dict(output=str(output), report=str(report), applied=len(plan),
                still_needs_review=[dict(name=r['name'], columns=[_one(c) for c in r['issue_columns']])
                                        for r in remaining])


def export_fasta(project_path, output_dir, allow_single=False, suffix=None):
    """One FASTA per marker from a saved project (stale contigs and, unless allowed, single reads are refused).

    Files are named ``<marker><suffix>.fasta``. By default the suffix is ``.reviewed`` when the project file
    is a reviewed copy (``*.reviewed.json``), so the batch run's unedited FASTA is never overwritten.
    """
    project_path = Path(project_path).expanduser()
    if suffix is None:
        suffix = '.reviewed' if project_path.stem.endswith('.reviewed') else ''
    project = load_project(project_path)
    rows = triage(project)
    skip = {r['key'] for r in rows if r['category'] == SINGLE_READ and not allow_single}
    keep = {k: c for k, c in project.contigs.items() if k not in skip}
    if not keep:
        raise ValueError('Nothing to export (single-read contigs need allow_single).')
    output_dir = Path(output_dir).expanduser()
    output_dir.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=output_dir) as tmp:
        names = export_markers(SimpleNamespace(contigs=keep), Path(tmp), allow_single=True)
        files = []
        for name in names:
            final = name[:-len('.fasta')] + suffix + '.fasta'
            shutil.move(str(Path(tmp) / name), str(output_dir / final))
            files.append(final)
    return dict(fasta=files, single_read_not_exported=sorted(r['name'] for r in rows if r['key'] in skip))
