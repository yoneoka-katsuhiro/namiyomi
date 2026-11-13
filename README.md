# namiyomi ver 0.2.0

A focused desktop editor for Sanger sequencing: inspect AB1 chromatograms, trim read ends, assemble selected reads, review the consensus, and export FASTA.

namiyomi keeps the sequence and its supporting traces in one workspace. It is designed for routine, small-scale Sanger work rather than as an integrated bioinformatics suite. Its features are **deliberately limited to what this work needs**, so there is little to configure or learn before you start. Assembly and editing run locally. Online AI review is optional. The application and source code are free under the MIT License.

## Download namiyomi here

- **[Download for macOS — Apple Silicon, macOS 15+](https://github.com/yoneoka-katsuhiro/namiyomi/releases/download/v0.2.0/namiyomi-0.2.0-macOS-arm64-development.zip)**
- **[Download for Windows — Windows 10/11 x64, experimental](https://github.com/yoneoka-katsuhiro/namiyomi/releases/download/v0.2.0/namiyomi-0.2.0-Windows-x64.zip)**

These ZIPs contain the ready-to-run application. No Python installation is required. The files in the repository are for developers; you do not need to download them to use namiyomi.

[Release notes](https://github.com/yoneoka-katsuhiro/namiyomi/releases/tag/v0.2.0) · [日本語](README-ja.md) · [Build and distribution guide](DISTRIBUTION.md)

![The consensus sequence above its aligned source chromatograms in namiyomi on macOS. The displayed data are synthetic.](docs/images/namiyomi-contig-macos.png)

*Review the consensus alongside its supporting chromatograms. Base calls, quality bars, and review flags share the same sequence coordinates. This screenshot uses the synthetic demonstration data supplied with the project.*

## Download and open

| Package | System | Release status |
| --- | --- | --- |
| `namiyomi-0.2.0-macOS-arm64-development.zip` | Apple Silicon, macOS 15 or later | Tested on the developer's Mac. Ad-hoc signed; not Developer ID signed or Apple notarized. Not an Intel or universal build. |
| `namiyomi-0.2.0-Windows-x64.zip` | Windows 10/11, x64 | Experimental. Full hands-on testing on a physical Windows PC has not been completed. Errors may occur. |

Download the appropriate ZIP from [Releases](https://github.com/yoneoka-katsuhiro/namiyomi/releases). Standalone packages include Python and their dependencies; a separate Python installation is not required.

**macOS:** extract the ZIP, open its `namiyomi` folder, then open `namiyomi.app`. You can move the application to an Applications folder if desired; keep the separate `demo` folder if you want to try the example project. This is an unnotarized development release. macOS may warn about or block the application; the package has not received Apple's notarization. See [Apple's explanation of Developer ID](https://developer.apple.com/developer-id/) before deciding whether to run it. Do not disable system security protections to install namiyomi.

**Windows:** extract the complete `namiyomi` folder, then open `namiyomi.exe`. Keep the `_internal` folder beside the executable; do not copy the executable out on its own. The executable is not Authenticode signed, so Windows may display a security warning. The Windows build workflow uses Windows Server 2022 for automated build and startup checks; consult [Actions](https://github.com/yoneoka-katsuhiro/namiyomi/actions) for the outcome of each run. Those checks do not replace hands-on testing of the complete workflow on Windows 10/11.

Version 0.2.0 is an early release. Back up your AB1 files and projects, and check exported sequences before using them in an analysis.

## What namiyomi does

- Imports AB1 files without assuming a particular filename pattern.
- Automatically trims low-quality ends on import, with adjustable settings and manual left/right cuts.
- Builds a contig from two or more selected reads and displays their aligned chromatograms beneath the consensus.
- Recalculates the linked consensus and its supporting-quality display when a source read is edited.
- Provides review flags, base editing, reverse complement, Undo/Redo, and sequence/FASTA copying.
- Finds consensus subsequences on either strand, with exact matching by default and an adjustable 80–100% identity threshold.
- Saves editable projects, writes local recovery snapshots, and exports FASTA without modifying the original AB1 files.

The interface deliberately omits cloning design, restriction mapping, BLAST, phylogenetics, and genome-scale assembly. Its scope is the work between opening Sanger reads and obtaining a reviewed sequence. The regular editor works offline without an account.

## Try the demonstration

Both application ZIPs include a `demo` folder with four synthetic AB1 files and a prepared `Demo.json` project. The same files are available in the repository's [`demo`](demo) directory. The sequences and trace signals were generated for demonstration; they are **not biological observations or instrument measurements**.

Keep the files in `demo` together, then choose **File > Import Previous Project** and open `demo/Demo.json` from the extracted application folder. You can inspect the prepared contigs, try Find, edit a base, and Undo the change. Use **Save Project As** to save your practice session separately. To practise assembly from scratch, start a new project and import the four demonstration AB1 files instead.

## A typical session

1. **Import the reads.** Choose **File > Import AB1** or drag AB1 files into the window. Reads remain separate in **AB1 files**; importing does not group them or build contigs. Automatic end trimming is enabled by default. Adjust **Settings > Trimming Settings** before importing if needed.
2. **Inspect the chromatograms.** Check trimmed ends and ambiguous calls. Gray regions are excluded from the retained sequence. Adjust **Base spacing** and **Trace height** at the bottom right. Use **Edit > Cut Left / Cut Right** for additional cuts.
3. **Select reads from one specimen and one locus.** Shift-click the relevant entries, then choose **Assembly > Build Contig (manual)**. Three or more reads are supported. The assembler checks both orientations; names and primer labels do not determine the alignment.
4. **Review the consensus.** Examine flagged positions against the traces. Click a base to edit it; Backspace removes selected calls. Source-read edits recalculate the linked consensus. Use **Review > Mark Reviewed** after checking a flagged call. A review acknowledgement does not improve its quality score.
5. **Save and export.** **Save Project** first asks for a filename, then updates that same file. **Export Project Copy** leaves the current save destination unchanged. Select a read or contig and choose **Export FASTA**, or use **Copy FASTA** from its right-click menu.

If an assembly needs to be redone, right-click its contig and choose **Return Reads to AB1**. The reads keep their edits and trimming. The removed contig's consensus-only edits and review flags are not transferred to them; Undo restores the whole operation.

### Optional filename-based assembly

**Assembly > Build Contig (AUTO)** groups unused reads by the part of the original filename before its first underscore. For example, `sample07_forward.ab1` and `sample07_reverse.ab1` form the group `sample07`. A successful group receives that name as its contig identifier.

Use AUTO only when each prefix identifies one specimen **and one locus**. Reads already used by an existing manual or AUTO contig are excluded, and existing contigs are not rebuilt. Use the manual command when filenames do not follow this convention. A single-read group can form a contig, but it does not provide independent multi-read confirmation.

### AI Review

Use **AI Review** to check flagged bases against their chromatograms.

1. **Start the review.** Open **AI Review > Review Start**. namiyomi uses the same filename grouping as AUTO to assemble unused reads. It then checks flagged positions on your device. Where the traces give a clear answer, it can fix a base or remove an extra base call. The default, **This device only**, needs no API key.
2. **Ask AI about the remaining flags (optional).** Before starting, open **AI Review > Settings**, choose OpenAI, Claude, or Gemini, and enter an API key. namiyomi applies a confident answer only if the trace also supports it. Manual edits and double peaks already recorded as IUPAC codes are preserved.

   The service receives the bases and quality scores near each position and a short excerpt of the signal as numbers. The original AB1 files, chromatogram images, and sample names stay on your device.
3. **See what needs a closer look.** In **AI Review Report**, **Flags** counts the positions still open. Click a contig name to jump to a flag in the chromatogram. Click the same name again to visit its next red flag.

![AI Review Settings with This device only selected.](docs/images/namiyomi-ai-settings-macos.jpg)

![AI Review Report for the included demo.](docs/images/namiyomi-ai-report-macos.jpg)

- **Red flag:** AI reviewed the position, but namiyomi left it open.
- **Black flag:** AUTO flagged the position. The mark stays after an automatic correction, and a corrected position no longer counts in **Flags**.
- **Teal underline:** namiyomi applied an automatic decision at that position.

Use **Edit > Undo** to undo the assembly and review together. See the [trace review guide](docs/TRACE_REVIEW.md) for how decisions are made.

### Process a folder of AB1 files from the command line

Use the batch command to trim reads, assemble them by sample, and export FASTA in one run. For example, `sample07_F.ab1` and `sample07_R.ab1` belong to the same group because the names match before the first `_`. Use each prefix for one specimen and one locus. If your filenames do not follow this pattern, supply a CSV to define the groups.

**Try the included demo** with Python 3.10 or later and the extracted source code ZIP. In a terminal, open the folder that contains `README.md` and run these two lines in order. If your system uses `python3` (Mac) or `py -3` (Windows), use that in place of `python`. [Download the v0.2.0 source code ZIP](https://github.com/yoneoka-katsuhiro/namiyomi/archive/refs/tags/v0.2.0.zip)

```sh
python -m pip install -r requirements/batch.txt
python src/batch.py --input demo --output demo-results
```

Install the dependency once. For your own data, replace `demo` with the folder containing your AB1 files and `demo-results` with a new output folder. Your original AB1 files stay intact.

Look in `demo-results` for:

- **`auto.fasta`** with the assembled sequences from the demo.
- **`project.json`**, which you can open in namiyomi with **File > Import Previous Project** to inspect the traces.
- **`summary.tsv`** with the assembly results and contigs that need review.

Add `--review` to check flagged positions against their traces on your device and save the decisions to `review.tsv`. These steps require no API key or MCP setup (the [MCP review server](#review-and-trim-from-an-mcp-client) is described below). For CSV grouping and optional online AI review, see the [batch guide](docs/BATCH.md#command-line).

### Review and trim from an MCP client

This is a different feature from **AI Review**, and you do not need it to use AI Review. The table shows how they differ.

| | **AI Review** (menu in the app) | **MCP review server** |
| --- | --- | --- |
| Who reviews | namiyomi's local rules, then optionally OpenAI, Claude or Gemini through **your API key** | An MCP client you already use, such as Claude Code or Claude Desktop, with no API key in namiyomi |
| How it runs | Click **AI Review > Review Start** in the app | You ask the client in chat, for example "run namiyomi on my reads, review the flags, then trim the ends" |
| What the AI model sees | Bases, quality scores and a short numeric signal excerpt near each flag | Sequences, file names and chromatogram images, through the tools below |
| Setup | None | `pip install -r requirements/mcp.txt` and register `src/mcp_review_server.py` with the client |
| Results | Written into the project and shown in the app | The same: written to a new project file that opens in the app, with flags the client left open shown in red |

Use AI Review for a quick check inside the app. Use the MCP server when you want to work with the reads in a conversation: batch a folder, look at individual traces, and trim contig ends by a rule you describe. Both settle the clear flags with the same local trace rules first.

1. **Install** (once): `python -m pip install -r requirements/mcp.txt`.
2. **Register the server** with your client. For Claude Code:
   ```sh
   claude mcp add namiyomi-review -- /path/to/python /path/to/namiyomi/src/mcp_review_server.py
   ```
   Other clients take the same command and argument in their MCP settings.
3. **Ask the client to work on a folder**, and give it an output folder that does not exist yet. It can run the batch, show flagged traces, and propose end trims. It changes nothing until you agree to a rule.
4. **Open the result.** Edits are saved to a new `project.reviewed.json` with a `review_report.md`; your original project and FASTA stay untouched. Open it with **File > Import Previous Project**. Columns are numbered from 1, as in the app.

The server makes no network calls itself, but the connected client sees the sequences and trace images it asks for, so register it only with a client you trust with your data. The smaller `src/mcp_server.py` is separate: it only reports anonymous quality counts and cannot edit or export anything. See the [batch guide](docs/BATCH.md#review-and-trim-from-an-mcp-client) for the full tool list.

### Find a sequence on either strand

![Orange highlighting marks a sequence match in a synthetic consensus, with its source chromatograms below.](docs/images/namiyomi-find-macos.png)

*Enter a DNA sequence in Find to highlight matching consensus sites. The search also checks the reverse complement. These demonstration traces are synthetic.*

The default is 100% identity. **Settings > Find Settings** permits 80–100% identity for substitution-only, ungapped matching. IUPAC letters match literally, not as wildcard sets. Search includes edited consensus calls and insertions, excludes deleted calls, and never changes the sequence. Use the arrow buttons or Enter / Shift+Enter to move between matches; Esc clears the query.

### Useful shortcuts

| Action | macOS | Windows |
| --- | --- | --- |
| Save project | ⌘S | Ctrl+S |
| Undo / Redo | ⌘Z / ⇧⌘Z | Ctrl+Z / Ctrl+Shift+Z |
| Build Contig (manual) | ⌘B | Ctrl+B |
| Copy sequence / Copy FASTA | ⌘C / ⇧⌘C | Ctrl+C / Ctrl+Shift+C |
| Find | ⌘F | Ctrl+F |
| Next / previous review flag | ⌘E / ⇧⌘E | Ctrl+E / Ctrl+Shift+E |

Right-click the sidebar for common actions. Double-click a sequence name to rename it. Double-click **Sequences** or **Length (bp)** to sort reads and contigs, then repeat to reverse the order. Delete, return-to-AB1, and export-with-flags confirmations offer **Don't ask again for this action**; **Settings > Confirmations** restores each prompt independently.

## Choosing a tool

namiyomi does one job: review AB1 chromatograms, build contigs, correct the consensus while looking at the traces, and export FASTA. The alternatives fall into four groups. The table shows what namiyomi and each group can and cannot do.

| Software | Operating systems | What it can do | What it cannot do |
| --- | --- | --- | --- |
| **namiyomi** | macOS (Apple Silicon, 15 or later)<br>Windows 10/11 (x64, experimental) | Sanger assembly and editing. Contigs built manually or automatically. A minimal, intuitive design that is easy to understand, optimized for research and teaching use. | Complex operations such as NGS read assembly and reference mapping are out of scope. |
| **Commercial Sanger suites** | **Sequencher:** Windows (the current version does not support 64-bit macOS)<br>**Geneious Prime:** Windows, macOS, Linux<br>**ChromasPro:** Windows only | Sanger assembly and editing, plus broader analyses (depending on the product: heterozygote and mutation detection, reference mapping, cloning tools, NGS), batch assembly with configurable naming rules (in several products), and vendor support. | They cannot be used free of charge for regular research work; a purchased or subscription license is required. You also have to find the steps you need among many general-purpose tools. |
| **Free chromatogram viewers** | **Chromas:** Windows only | Viewing, editing, trimming, searching, and exporting a single chromatogram. | It cannot assemble reads into a contig; its developer directs users to the commercial ChromasPro for that. |
| **Free bioinformatics toolkits** | **UGENE:** Windows, macOS, Linux | Many analysis tools in one program, including a Sanger reads editor and CAP3 assembly. | It is not dedicated to Sanger review; you choose the steps you need from a large set of general-purpose tools. |
| **R and command-line tools** | **sangeranalyseR:** R on Windows, macOS, Linux<br>**Tracy:** macOS, Linux (Bioconda, Homebrew, containers) | Scripted, reproducible batch processing and reports. | They cannot be used without R or the command line, and you cannot correct individual bases while looking at the traces (sangeranalyseR's documentation directs manual editing to another tool). |

Choose one of the other tools when you need what it adds. namiyomi has not been benchmarked against them for speed or accuracy. The [detailed comparison](docs/COMPARISON.md) lists each application's documented features, platforms, and licensing, with links to the developers' documentation.

## Important limits and data handling

- **Assembly requires related, overlapping reads.** The greedy overlap assembler is for small Sanger contigs, not genome assembly or repeat resolution. Defaults are a minimum overlap of 100 comparable bases and 90% identity. A filename or overlap cannot establish specimen identity, locus identity, or absence of contamination.
- **Review flags identify positions to inspect.** Conflicting high-quality calls and internal alignment gaps can remain `N`. **Support Q** is the highest usable original Phred score supporting an automatic call, not a calibrated combined confidence probability. Manually edited source calls do not acquire instrument Phred scores.
- **Relative orientation is not biological annotation.** Orientation unification uses earlier sufficiently similar contigs as references. The first contig in each matching locus sets the relative direction; namiyomi does not determine a gene's biological forward strand.
- **A project is not a self-contained AB1 archive.** Project JSON stores source-file paths and checksums, not the chromatograms themselves. Keep the original AB1 files available. When moving a project, include the referenced AB1 files and preserve their relative folder structure. Changed or missing sources can prevent reopening.
- **Recovery is a safeguard, not a backup.** Unsaved edits receive local recovery snapshots approximately every 30 seconds. A crash can still lose the most recent edits and Undo history. Save regularly and keep independent backups.

## License, citation, and feedback

namiyomi is released under the [MIT License](LICENSE.txt). Bundled dependencies retain their own licenses and are accompanied by third-party notices.

Suggested citation:

> Yoneoka, K. (2026). *Namiyomi: Sanger chromatogram review and contig assembly* (Version 0.2.0) [Computer software]. GitHub. https://github.com/yoneoka-katsuhiro/namiyomi

In text, cite **Yoneoka (2026)** or **(Yoneoka, 2026)**.

Machine-readable metadata is available in [CITATION.cff](CITATION.cff). Record the software version, processing settings, and manual review decisions when describing your workflow. No DOI has been assigned.

Report a problem through [GitHub Issues](https://github.com/yoneoka-katsuhiro/namiyomi/issues). Include the operating system, namiyomi version, steps to reproduce, and the error message or a screenshot. Do not upload confidential or unpublished sequence data unless you intend to make it public.

Developed by Katsuhiro Yoneoka, Mishima, Shizuoka, Japan.
