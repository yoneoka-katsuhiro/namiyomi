import tempfile
import unittest
from pathlib import Path

from scripts.create_demo import create_demo

try:
    import mcp  # noqa: F401
    HAVE_MCP = True
except ImportError:  # the optional batch/MCP requirements are not installed
    HAVE_MCP = False
if HAVE_MCP:
    import mcp_review_server  # a broken server must fail the tests, not be skipped


@unittest.skipUnless(HAVE_MCP, 'mcp package not installed (requirements/batch.txt)')
class McpToolTests(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        root = Path(self._tmp.name)
        self.demo, self.out = root / 'demo', root / 'out'
        create_demo(self.demo)

    def tearDown(self):
        self._tmp.cleanup()

    def test_registered_tools(self):
        import asyncio
        names = {t.name for t in asyncio.run(mcp_review_server.mcp.list_tools())}
        self.assertEqual(names, {'scan_folder', 'run_batch', 'list_contigs', 'get_contig_details',
                                 'get_column_qualities', 'describe_contig_ends', 'suggest_quality_trim',
                                 'find_sequence', 'render_trace', 'apply_review_edits', 'export_fasta'})

    def test_trim_prompt(self):
        import asyncio
        prompts = {p.name for p in asyncio.run(mcp_review_server.mcp.list_prompts())}
        self.assertEqual(prompts, {'trim_contig_ends'})
        text = mcp_review_server.trim_contig_ends('/x/project.json', 'cut before 4 Ns')
        self.assertIn('describe_contig_ends', text)
        self.assertIn('"cut before 4 Ns"', text)
        self.assertIn('Do not apply anything', text)

    def test_human_review_is_default_and_makes_no_edits(self):
        self.assertEqual(mcp_review_server.scan_folder(str(self.demo))['ab1_files'], 4)
        summary = mcp_review_server.run_batch(str(self.demo), str(self.out))
        self.assertEqual(summary['review_mode'], 'human')
        self.assertIn('do not edit anything unless', summary['next_steps'])
        self.assertFalse((self.out / 'project.reviewed.json').exists())
        with self.assertRaises(ValueError):
            mcp_review_server.run_batch(str(self.demo), str(self.out), review='robot')
        with self.assertRaisesRegex(ValueError, 'already exists'):   # names the folder, suggests another
            mcp_review_server.run_batch(str(self.demo), str(self.out))

    def test_claude_review_flow(self):
        summary = mcp_review_server.run_batch(str(self.demo), str(self.out), review='claude')
        self.assertIn('render_trace', summary['next_steps'])
        project = str(self.out / 'project.json')
        flagged = mcp_review_server.list_contigs(project, 'needs_review')[0]
        details = mcp_review_server.get_contig_details(project, flagged['name'])
        col = details['flagged'][0]['column']
        self.assertEqual(mcp_review_server.render_trace(project, flagged['name'], col).data[:4], b'\x89PNG')
        result = mcp_review_server.apply_review_edits(project, str(self.out / 'project.reviewed.json'),
            [dict(contig=flagged['name'], column=col, action='set', value='C', reason='clean C peak')])
        self.assertEqual(result['still_needs_review'], [])
        files = mcp_review_server.export_fasta(str(self.out / 'project.reviewed.json'), str(self.out / 'final'))
        self.assertEqual(files['fasta'], ['auto.reviewed.fasta'])


if __name__ == '__main__':
    unittest.main()
