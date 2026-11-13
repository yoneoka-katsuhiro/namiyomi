import tempfile
import unittest
from pathlib import Path

import batch
import review
from scripts.create_demo import create_demo

PNG = b'\x89PNG\r\n\x1a\n'
try:
    import matplotlib  # noqa: F401
    HAVE_MPL = True
except ImportError:
    HAVE_MPL = False


class ReviewTests(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        root = Path(self._tmp.name)
        demo = root / 'demo'
        create_demo(demo)
        self.out = root / 'out'
        self.summary = batch.run_batch(demo, self.out)
        self.project_path = self.out / 'project.json'
        self.reviewed = self.out / 'project.reviewed.json'
        self.flag = self.summary['contigs'][0]['issue_columns'][0]  # Demo01 has one conflict (0-based)
        self.col = self.flag + 1                                     # the same column as the tools number it

    def tearDown(self):
        self._tmp.cleanup()

    def test_list_and_details(self):
        project = review.load_project(self.project_path)
        self.assertEqual([r['name'] for r in review.list_contigs(project, 'needs_review')], ['Demo01'])
        self.assertEqual([r['name'] for r in review.list_contigs(project, 'clean')], ['Demo02'])
        with self.assertRaises(ValueError):
            review.list_contigs(project, 'bogus')
        d = review.contig_details(project, 'Demo01')
        self.assertEqual(len(d['sources']), 2)
        flag = d['flagged'][0]
        self.assertEqual(flag['column'], self.col)  # columns are numbered from 1, as in the app
        self.assertEqual(flag['reason'], 'conflict')
        self.assertEqual(len(flag['reads']), 2)
        self.assertTrue(all(r['base'] and r['phred'] is not None for r in flag['reads']))
        for r in flag['reads']:
            heights = r['peak_heights']
            self.assertEqual(sorted(heights), list('ACGT'))
            self.assertEqual(max(heights, key=heights.get), r['base'])  # called base has the tallest peak
        self.assertEqual(review.contig_details(project, d['key'])['name'], 'Demo01')  # key works too
        with self.assertRaises(ValueError):
            review.contig_details(project, 'nope')

    def test_find_sequence(self):
        project = review.load_project(self.project_path)
        seq = next(c for c in project.contigs.values() if c.sequence_name == 'Demo02').sequence()
        hits = review.find_sequence(project, 'Demo02', seq[100:120])
        self.assertEqual(hits[0]['start'], 101)  # positions are numbered from 1
        self.assertEqual(hits[0]['strand'], '+')

    @unittest.skipUnless(HAVE_MPL, 'matplotlib not installed')
    def test_render_trace(self):
        project = review.load_project(self.project_path)
        for ref, col in (('Demo01', self.col), ('Demo01', 1), ('Demo02', 300)):
            png = review.render_trace(project, ref, col)
            self.assertTrue(png.startswith(PNG) and len(png) > 5000)
        with self.assertRaises(ValueError):
            review.render_trace(project, 'Demo01', 10 ** 6)
        with self.assertRaises(ValueError):
            review.render_trace(project, 'Demo01', 0)  # numbering starts at 1

    def test_apply_edits_writes_copy_and_report(self):
        before = self.project_path.read_bytes()
        result = review.apply_review_edits(self.project_path, self.reviewed, [
            dict(contig='Demo01', column=self.col, action='set', value='a', reason='Clean A peak in F')])
        self.assertEqual(self.project_path.read_bytes(), before)  # input untouched
        self.assertEqual(result['still_needs_review'], [])
        project = review.load_project(self.reviewed)
        key, contig = review.resolve_contig(project, 'Demo01')
        orig = review.resolve_contig(review.load_project(self.project_path), 'Demo01')[1]
        self.assertEqual(contig.call(self.flag), 'A')
        self.assertEqual(contig.issues(), [])
        other = [i for i in range(len(contig.columns)) if i != self.flag]
        self.assertTrue(all(contig.call(i) == orig.call(i) for i in other))
        report = Path(result['report']).read_text(encoding='utf-8')
        self.assertIn('Clean A peak in F', report)
        self.assertIn('Left unresolved', report)
        review.apply_review_edits(self.project_path, self.out / 'second.json', [
            dict(contig='Demo01', column=self.col, action='confirm', reason='kept N')],
            report_path=result['report'])
        self.assertEqual(Path(result['report']).read_text(encoding='utf-8').count('## Claude'), 2)

    def test_confirm_marks_reviewed_without_changing_call(self):
        review.apply_review_edits(self.project_path, self.reviewed, [
            dict(contig='Demo01', column=self.col, action='confirm', reason='N is correct here')])
        contig = review.resolve_contig(review.load_project(self.reviewed), 'Demo01')[1]
        self.assertEqual(contig.call(self.flag), 'N')
        self.assertNotIn(self.flag, contig.issues())
        self.assertEqual(contig.black_marks[str(self.flag)], False)

    def test_decisions_use_the_in_app_review_store(self):
        review.apply_review_edits(self.project_path, self.reviewed, [
            dict(contig='Demo01', column=self.col, action='set', value='c', reason='clean C peak\nin F')])
        contig = review.resolve_contig(review.load_project(self.reviewed), 'Demo01')[1]
        key = str(self.flag)
        self.assertEqual(contig.ai_review[key], 'C')
        self.assertNotIn(key, contig.overrides)          # hand edits stay the user's
        self.assertTrue(contig.ai_decided(self.flag))
        self.assertTrue(contig.ai_notes[key].startswith('Claude (MCP): clean C peak'))
        self.assertNotIn('\n', contig.ai_notes[key])
        self.assertEqual(contig.call(self.flag), 'C')

    def test_remove_column_and_long_reason(self):
        review.apply_review_edits(self.project_path, self.reviewed, [
            dict(contig='Demo01', column=self.col, action='set', value='-', reason='x' * 500)])
        contig = review.resolve_contig(review.load_project(self.reviewed), 'Demo01')[1]
        self.assertEqual(contig.ai_review[str(self.flag)], '')
        self.assertEqual(contig.call(self.flag), '')
        self.assertLessEqual(len(contig.ai_notes[str(self.flag)]), review.NOTE_LIMIT)

    def test_unresolved_shows_as_red_flag_in_app(self):
        review.apply_review_edits(self.project_path, self.reviewed, [], unresolved=[
            dict(contig='Demo01', column=self.col, reason='overlapping peaks')])
        contig = review.resolve_contig(review.load_project(self.reviewed), 'Demo01')[1]
        self.assertTrue(contig.ai_unresolved(self.flag))
        self.assertEqual(contig.flag_color(self.flag), 'red')
        self.assertIn('overlapping peaks', contig.ai_notes[str(self.flag)])

    def test_hand_edited_columns_are_not_overwritten(self):
        total = len(review.resolve_contig(review.load_project(self.project_path), 'Demo01')[1].columns)
        review.apply_review_edits(self.project_path, self.reviewed, [
            dict(contig='Demo01', keep_from=2, keep_to=total, action='trim', reason='test trim')])
        with self.assertRaisesRegex(ValueError, 'edited by hand'):   # column 1 was trimmed
            review.apply_review_edits(self.reviewed, self.out / 'again.json', [
                dict(contig='Demo01', column=1, action='set', value='A', reason='no')])

    def test_columns_are_numbered_from_one(self):
        for bad in (0, -1, 'x'):
            with self.assertRaises(ValueError, msg=str(bad)):
                review.apply_review_edits(self.project_path, self.reviewed, [
                    dict(contig='Demo01', column=bad, action='set', value='A', reason='x')])
        d = review.contig_details(review.load_project(self.project_path), 'Demo01')
        self.assertEqual(d['flagged'][0]['column'], self.flag + 1)
        result = review.apply_review_edits(self.project_path, self.reviewed, [], unresolved=[
            dict(contig='Demo01', column=self.col, reason='r')])
        self.assertEqual(result['still_needs_review'][0]['columns'], [self.col])
        contig = review.resolve_contig(review.load_project(self.reviewed), 'Demo01')[1]
        self.assertTrue(contig.ai_unresolved(self.flag))     # stored 0-based, like the app's own review

    def test_invalid_edits_change_nothing(self):
        good = dict(contig='Demo01', column=self.col, action='set', value='A', reason='ok')
        bad_cases = [dict(good, value='Z'), dict(good, value='AC'), dict(good, value='R'), dict(good, reason=' '),
                     dict(good, column=10 ** 6), dict(good, column=0), dict(good, action='delete'), dict(good, contig='nope')]
        for bad in bad_cases:
            with self.assertRaises(ValueError, msg=str(bad)):
                review.apply_review_edits(self.project_path, self.reviewed, [good, bad])
            self.assertFalse(self.reviewed.exists())
        with self.assertRaises(ValueError):
            review.apply_review_edits(self.project_path, self.project_path, [good])
        with self.assertRaises(ValueError):
            review.apply_review_edits(self.project_path, self.reviewed, [])

    def test_unresolved_flags_are_listed_in_report(self):
        result = review.apply_review_edits(self.project_path, self.reviewed, [], unresolved=[
            dict(contig='Demo01', column=self.col, reason='overlapping peaks')])
        self.assertEqual(result['applied'], 0)
        self.assertEqual(result['still_needs_review'], [dict(name='Demo01', columns=[self.col])])
        report = Path(result['report']).read_text(encoding='utf-8')
        self.assertIn(f'**Demo01** column {self.col}: overlapping peaks', report)
        with self.assertRaises(ValueError):
            review.apply_review_edits(self.project_path, self.out / 'x.json', [],
                                      unresolved=[dict(contig='nope', column=1)])

    def test_export_reviewed_does_not_overwrite_batch_fasta(self):
        before = (self.out / 'auto.fasta').read_text()
        review.apply_review_edits(self.project_path, self.reviewed, [
            dict(contig='Demo01', column=self.col, action='set', value='A', reason='clean A')])
        result = review.export_fasta(self.reviewed, self.out)
        self.assertEqual(result['fasta'], ['auto.reviewed.fasta'])
        self.assertEqual((self.out / 'auto.fasta').read_text(), before)
        self.assertNotEqual((self.out / 'auto.reviewed.fasta').read_text(), before)
        self.assertEqual(sorted(p.name for p in self.out.iterdir() if p.is_dir()), [])  # temp dir removed

    def test_column_qualities(self):
        project = review.load_project(self.project_path)
        q = review.column_qualities(project, 'Demo01')
        self.assertEqual(len(q['rows']), q['columns'])
        row = next(r for r in q['rows'] if r['column'] == self.col)
        self.assertEqual(q['rows'][0]['column'], 1)
        self.assertEqual(len(row['reads']), 2)
        self.assertEqual(len(review.column_qualities(project, 'Demo01', 10, 19)['rows']), 10)
        with self.assertRaises(ValueError):
            review.column_qualities(project, 'Demo01', 20, 10)
        with self.assertRaises(ValueError):
            review.column_qualities(project, 'Demo01', 0, 10)

    def test_trim_deletes_outside_columns(self):
        orig = review.resolve_contig(review.load_project(self.project_path), 'Demo01')[1]
        last = len(orig.columns) - 1
        result = review.apply_review_edits(self.project_path, self.reviewed, [
            dict(contig='Demo01', action='trim', keep_from=6, keep_to=last - 9, reason='rolling mean Q < 20')])
        contig = review.resolve_contig(review.load_project(self.reviewed), 'Demo01')[1]
        expected = ''.join(orig.call(i) for i in range(5, last - 9)).replace('-', '')
        self.assertEqual(contig.sequence(), expected)
        self.assertIn('trimmed to columns 6-', Path(result['report']).read_text(encoding='utf-8'))
        for bad in (dict(keep_from=1, keep_to=last + 1), dict(keep_from=10, keep_to=5), dict(keep_from=0, keep_to=5),
                    dict(keep_from=1, keep_to=last + 2)):
            with self.assertRaises(ValueError, msg=str(bad)):
                review.apply_review_edits(self.project_path, self.out / 'bad.json',
                                          [dict(contig='Demo01', action='trim', reason='x', **bad)])

    def test_rolling_keep(self):
        self.assertEqual(review._rolling_keep([5, 5, 30, 30, 30], 2, 20, 5), 2)
        self.assertEqual(review._rolling_keep([5, None, 30, 30], 2, 20, 4), 1)  # None left out of the mean
        self.assertEqual(review._rolling_keep([5, 5, 5, 30], 1, 20, 2), 2)      # never past the limit
        self.assertEqual(review._rolling_keep([30, 5], 1, 20, 2), 0)

    def test_suggest_quality_trim(self):
        project = review.load_project(self.project_path)
        lenient = review.suggest_quality_trim(project, window=10, min_mean_phred=0)
        self.assertEqual({r['name'] for r in lenient}, {'Demo01', 'Demo02'})
        self.assertTrue(all(r['edit'] is None for r in lenient))
        strict = review.suggest_quality_trim(project, 'Demo01', window=5, min_mean_phred=99,
                                             single_read_ends_only=True)[0]
        self.assertEqual([strict['keep_from'], strict['keep_to']], strict['overlap'])  # stops at the overlap
        free = review.suggest_quality_trim(project, 'Demo01', window=5, min_mean_phred=99,
                                           single_read_ends_only=False)[0]
        self.assertIn('error', free)
        before = review.resolve_contig(project, 'Demo01')[1].sequence()
        self.assertEqual(review.resolve_contig(review.load_project(self.project_path), 'Demo01')[1].sequence(),
                         before)  # dry run changed nothing
        if strict['edit']:
            result = review.apply_review_edits(self.project_path, self.reviewed, [strict['edit']])
            self.assertEqual(result['applied'], 1)
        with self.assertRaises(ValueError):
            review.suggest_quality_trim(project, window=0)

    def test_n_run_cut(self):
        calls = list('AANCNNNAA') + ['', 'N']
        self.assertEqual(review._n_run_cut(calls, 3, range(len(calls))), 4)
        self.assertEqual(review._n_run_cut(calls, 1, range(len(calls))), 2)
        self.assertIsNone(review._n_run_cut(calls, 4, range(len(calls))))
        self.assertEqual(review._n_run_cut(list('ANN') + ['-'] + list('NA'), 3, range(6)), 1)  # gaps skipped

    def _project_with_n_runs(self):
        """Demo02 with an NNNN run in its right single-read end and an NN pair in the left one."""
        project = review.load_project(self.project_path)
        _, contig = review.resolve_contig(project, 'Demo02')
        a, b = review._Ends(contig).overlap
        for col in range(b + 20, b + 24):
            contig.overrides[str(col)] = 'N'
        for col in (a - 30, a - 29):
            contig.overrides[str(col)] = 'N'
        return project, contig, a, b

    def test_suggest_trim_max_n_run(self):
        project, contig, a, b = self._project_with_n_runs()
        r = review.suggest_quality_trim(project, 'Demo02', min_mean_phred=None, max_n_run=4)[0]
        self.assertEqual((r['keep_from'], r['keep_to']), (1, b + 20))
        self.assertIn('run of 4 or more Ns', r['edit']['reason'])
        r = review.suggest_quality_trim(project, 'Demo02', min_mean_phred=None, max_n_run=2)[0]
        self.assertEqual((r['keep_from'], r['keep_to']), (a - 27, b + 20))
        self.assertEqual(r['n_count'], 0)
        both = review.suggest_quality_trim(project, 'Demo02', window=5, min_mean_phred=99, max_n_run=4)[0]
        self.assertEqual((both['keep_from'], both['keep_to']), (a + 1, b + 1))  # stricter rule wins at each end
        with self.assertRaises(ValueError):
            review.suggest_quality_trim(project, min_mean_phred=None)
        with self.assertRaises(ValueError):
            review.suggest_quality_trim(project, max_n_run=0)

    def test_describe_contig_ends(self):
        project, contig, a, b = self._project_with_n_runs()
        d = review.describe_contig_ends(project, 'Demo02')[0]
        self.assertEqual(d['overlap'], [a + 1, b + 1])
        self.assertEqual(len(d['reads']), 2)
        self.assertEqual(d['right_end']['longest_n_run'], 4)
        self.assertEqual(d['left_end']['n_count'], 2)
        opts = {o['label']: o for o in d['options']}
        self.assertEqual(opts['keep everything']['length'], len(contig.sequence()))
        self.assertEqual(opts['cut before a run of >=4 N']['keep_to'], b + 20)
        self.assertEqual((opts['region covered by two reads']['keep_from'],
                          opts['region covered by two reads']['keep_to']), (a + 1, b + 1))
        custom = review.describe_contig_ends(project, 'Demo02', [dict(label='x', min_mean_phred=None, max_n_run=2),
                                                                 dict(label='bad', max_n_run=0)])[0]['options']
        self.assertEqual(custom[0]['keep_from'], a - 27)
        self.assertIn('error', custom[1])
        self.assertEqual(len(review.describe_contig_ends(project)), 2)

    def test_export_fasta_from_project(self):
        out = self.out / 'again'
        result = review.export_fasta(self.project_path, out)
        self.assertEqual(result['fasta'], ['auto.fasta'])
        self.assertEqual((out / 'auto.fasta').read_text(), (self.out / 'auto.fasta').read_text())


if __name__ == '__main__':
    unittest.main()
