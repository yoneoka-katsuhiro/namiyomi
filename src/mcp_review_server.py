"""Optional full-access stdio MCP server: batch tool, trace review and end trimming for an MCP client.

This is the opt-in companion to ``mcp_server.py`` (the read-only, anonymous quality summary).
Unlike that server it returns sequences, read names and trace images to the connected host, and it can
write new project files, so register it only with a client you trust with the data. The server itself
makes no network calls. Run with ``python src/mcp_review_server.py`` (see docs/BATCH.md).
Needs ``pip install -r requirements/mcp.txt``.
"""
from pathlib import Path

from mcp.server.fastmcp import FastMCP, Image

import batch
import review

mcp = FastMCP('namiyomi')

HUMAN_NEXT = ('Review mode is "human": by default do not edit anything unless the user explicitly asks you to '
              'review or clean the reads. Report the triage to the user and give them '
              'the project path to open in namiyomi (File > Import Previous Project), listing the needs_review, '
              'single_read and failed contigs for them to check against the traces.')
CLAUDE_NEXT = ('Review mode is "claude": for each needs_review contig call get_contig_details, then render_trace '
               'for every flagged column and judge the peaks. Apply only confident changes with '
               'apply_review_edits (writes a NEW reviewed project plus review_report.md); leave anything '
               'uncertain unedited, pass it as unresolved with a short reason so it is listed in the report, '
               'and list it for the user. Then export_fasta from the reviewed project (files are named '
               '<marker>.reviewed.fasta). Your reading of traces is a second opinion; say so.')
TRIM_NOTE = (' If the user asks about trimming contig ends, start with describe_contig_ends, agree on a rule, '
             'preview it with suggest_quality_trim and apply it only after they confirm.')


def _path(value):
    return Path(value).expanduser()


@mcp.tool()
def scan_folder(input_dir: str, sheet: str | None = None) -> dict:
    """Dry run: count the AB1 files in a folder and show the contig groups that would be built.

    Without ``sheet``, files are grouped by the text before the first underscore (like the app's
    Build Contig AUTO). With ``sheet`` (CSV with columns file,sample,marker) groups come from the sheet;
    files not in the sheet and sheet rows with no file are reported. Nothing is written.
    """
    return batch.run_batch(_path(input_dir), None, sheet=_path(sheet) if sheet else None, dry_run=True)


@mcp.tool()
def run_batch(input_dir: str, output_dir: str, sheet: str | None = None, allow_single: bool = False,
              review: str = 'human', min_quality: int | None = None, min_overlap: int | None = None,
              min_identity: float | None = None, local_rules: bool = True) -> dict:
    """Import, trim, group, assemble and triage every AB1 in input_dir; write results to output_dir.

    output_dir must be a different folder (originals are never modified). Writes project.json (open it
    in the app), summary.json/tsv and one FASTA per marker. Single-read contigs are exported only with
    allow_single. ``review`` is "human" (default: stop after triage, the user checks flags in the app) or
    "claude" (Claude reviews flags with render_trace and apply_review_edits). ``local_rules`` (default on)
    first settles the flags the app's own trace rules can decide, as the in-app review does; only the
    rest is left for you. No online model is called. Follow ``next_steps``.
    """
    if review not in ('human', 'claude'):
        raise ValueError('review must be "human" or "claude".')
    settings = {k: v for k, v in (('min_quality', min_quality), ('min_overlap', min_overlap),
                                  ('min_identity', min_identity)) if v is not None}
    out = _path(output_dir)
    if out.is_symlink() or out.exists():
        raise ValueError(f'{out} already exists. Choose a folder that does not exist yet (for example '
                         f'{out / "run1"}); existing results are never overwritten.')
    summary = batch.run_batch(_path(input_dir), out, sheet=_path(sheet) if sheet else None,
                              settings=settings, allow_single=allow_single,
                              review=local_rules)
    summary['review_mode'] = review
    summary['next_steps'] = (HUMAN_NEXT if review == 'human' else CLAUDE_NEXT) + TRIM_NOTE
    return summary


@mcp.tool()
def list_contigs(project_path: str, filter: str | None = None) -> list[dict]:
    """List contigs in a saved project.json, optionally filtered: clean, needs_review or single_read."""
    return review.list_contigs(review.load_project(_path(project_path)), filter)


@mcp.tool()
def get_contig_details(project_path: str, contig: str) -> dict:
    """Flagged columns of one contig (alignment columns, numbered from 1 like the app), with each source read's base call, Phred
    and peak_heights (raw A/C/G/T channel heights at that read's peak, already complemented for
    reverse reads so they match the contig orientation). None means the read does not cover the column.
    """
    return review.contig_details(review.load_project(_path(project_path)), contig)


@mcp.tool()
def get_column_qualities(project_path: str, contig: str, start: int | None = None,
                         end: int | None = None) -> dict:
    """Every alignment column in start..end (numbered from 1, inclusive; default: whole contig) with the consensus
    call and each read's [base, Phred]. Use it to pick a trim point when the user asks for quality-based
    end trimming (e.g. a rolling mean of Phred falling below a threshold).
    """
    return review.column_qualities(review.load_project(_path(project_path)), contig, start, end)


@mcp.tool()
def describe_contig_ends(project_path: str, contig: str | None = None,
                         options: list[dict] | None = None) -> list[dict]:
    """Start here for end trimming. For each contig (all unless ``contig`` is given): read layout, the
    region two reads cover, each single-read end (length, mean Phred, N count, longest N run, N columns),
    and what candidate rules would keep (range, length, N count). Nothing is changed.

    ``options`` overrides the compared rules: a list of {"label": str, plus any suggest_quality_trim
    argument (window, min_mean_phred, max_n_run, single_read_ends_only, keep_isolated_n), or
    "overlap_only": true}. The default compares keep-all, rolling means, N-run cuts and overlap-only.
    """
    return review.describe_contig_ends(review.load_project(_path(project_path)), contig, options)


@mcp.tool()
def suggest_quality_trim(project_path: str, contig: str | None = None, window: int = 10,
                         min_mean_phred: float | None = 20, single_read_ends_only: bool = True,
                         keep_isolated_n: bool = True, max_n_run: int | None = None) -> list[dict]:
    """Dry run of one end-trimming rule (all contigs unless ``contig`` is given); nothing is changed.

    Rolling mean: column quality = highest Phred among covering reads; from each end, columns are trimmed
    while the mean over the next ``window`` columns is below ``min_mean_phred`` (null turns this off).
    keep_isolated_n leaves single Ns out of the mean. N runs: ``max_n_run`` cuts each end just before
    the first run of that many consecutive Ns, walking outward from the two-read region. With both, the
    stricter cut wins. single_read_ends_only never cuts into the region two reads cover. Show the user
    the result, then pass the non-null ``edit`` values to apply_review_edits once they agree.
    Prefer this over computing a trim rule yourself from get_column_qualities.
    """
    return review.suggest_quality_trim(review.load_project(_path(project_path)), contig, window,
                                       min_mean_phred, single_read_ends_only, keep_isolated_n, max_n_run)


@mcp.prompt()
def trim_contig_ends(project_path: str, rule: str = '') -> str:
    """Guided end trimming of a batch result: compare rules, agree on one, apply it, export."""
    ask = f'The user\'s rule: "{rule}".' if rule else 'The user has not given a rule yet.'
    return f"""Help me trim the ends of the contigs in {project_path} with the namiyomi tools. {ask}

1. Call describe_contig_ends on the project. Summarise for me, in a short table per contig that has
   single-read ends: the two-read region, each single-read end's length, mean Phred, N count and longest
   N run, and what each compared rule would keep (length, N count). Leave out contigs where no rule
   changes anything, but say how many there were.
2. If I gave a rule, translate it into suggest_quality_trim arguments and say exactly how you read it
   (e.g. which ends, how isolated Ns are handled). If my rule cannot be expressed with those arguments,
   say so instead of approximating it silently. If I gave no rule, start conservative: recommend
   overlap-only (keeping just the region two reads cover) when it exists, but also lay out the
   alternatives from the comparison (rolling-mean Phred with a stricter threshold or longer window,
   cutting at runs of Ns, keeping everything) and explain the trade-off (bases kept versus Ns kept), so
   I can ask for a follow-up run. A rolling-mean rule scores each column by its best covering read and
   ignores isolated Ns, so on single-read tails with low Phred (below about 20) the default settings can
   change nothing; if the comparison shows that, say so and show a stricter variant.
3. Run suggest_quality_trim with the agreed arguments and show the proposed cuts. Do not apply anything
   until I confirm.
4. After I confirm, call apply_review_edits with the non-null edit values, writing project.reviewed.json
   next to the project (if that file already exists, use it as the input so earlier review edits are
   kept, and write to a new file name). Then export_fasta from the reviewed project and report the
   files written and the review_report.md entry.

Trimming is my decision; you are helping me apply it consistently. Never trim into the region two reads
cover unless I ask for that explicitly."""


@mcp.tool()
def find_sequence(project_path: str, contig: str, query: str, identity: int = 100) -> list[dict]:
    """Search a contig's consensus (both strands) for a DNA/IUPAC query; identity is 80-100 (%)."""
    return review.find_sequence(review.load_project(_path(project_path)), contig, query, identity)


@mcp.tool()
def render_trace(project_path: str, contig: str, column: int, window: int = 15) -> Image:
    """Image of every source read's chromatogram around one alignment column (highlighted).

    Base calls and Phred scores are shown above the peaks; reverse reads are shown complemented so
    all reads line up. Use it to judge a flagged column. ``window`` = columns shown on each side.
    """
    return Image(data=review.render_trace(review.load_project(_path(project_path)), contig, column, window),
                 format='png')


@mcp.tool()
def apply_review_edits(project_path: str, output_path: str, edits: list[dict],
                       unresolved: list[dict] | None = None) -> dict:
    """Apply review verdicts to a NEW copy of the project (the input is never overwritten).

    Each edit: {"contig": name, "column": int, "action": "set" or "confirm", "value": "A" (for set),
    "reason": "why, citing the peaks"}. Columns are numbered from 1, as in the app and review.tsv. "set"
    records a decision for the column (A, C, G, T, or "-" to remove it); "confirm" marks the flag reviewed
    without changing it. To trim a contig's ends: {"contig": name, "action": "trim",
    "keep_from": int, "keep_to": int, "reason": "the rule used"}; columns outside the range are deleted
    from the consensus. Trim only when the user asked for it, using their rule (suggest_quality_trim
    gives ready edits for rolling-mean Phred rules). Every edit is logged with its reason in review_report.md next to
    output_path. Do not edit columns you are unsure about; leave them for the user and pass them in
    ``unresolved`` as {"contig": name, "column": int, "reason": "why"}. Every flag still open is
    listed under "Left unresolved" in the report.
    """
    return review.apply_review_edits(_path(project_path), _path(output_path), edits, unresolved=unresolved)


@mcp.tool()
def export_fasta(project_path: str, output_dir: str, allow_single: bool = False) -> dict:
    """Write one FASTA per marker from a saved project (use the reviewed project after review).

    From a reviewed project (*.reviewed.json) the files are named <marker>.reviewed.fasta, so the
    batch run's unedited <marker>.fasta in the same folder is kept.
    """
    return review.export_fasta(_path(project_path), _path(output_dir), allow_single)


if __name__ == '__main__':
    mcp.run()
