# MitoPilot command-line mode: design

Date: 2026-11-25. Branch: `cli-mode` (off `devel`).

## Intent

Some users want a "productized" MitoPilot: less interaction, no Shiny app
for routine work, and no R console. Give them a `mitopilot` shell command
that runs a project from raw reads to GenBank files, fully scriptable, while
keeping a human in the loop for samples that need judgement.

## Goals

- One `mitopilot` command with subcommands, usable from a shell script.
- Clean samples go from reads to export with no human action.
- Samples that fail a pass rule are parked with a stated reason.
- Parked samples are cleared through a terminal walk-through that opens the
  app focused on the sample when eyes are needed, or takes a terminal-only
  decision when they are not (e.g. a failed assembly).
- Every setting editable without the app: guided editor, one-line commands,
  JSON file.
- App and command line share the same project database and the same rules.

## Non-goals

- Re-implementing deep edits in the terminal: gene-end adjustment, consensus
  building, scaffold join layout, alignment fixes. These stay in the app.
- Terminal plots or a full-screen terminal interface.
- Running `mitopilot` itself from inside the container image with no R on
  the host (possible later; not now).
- A user-defined rule language for auto-accept.

## Anti-goals

- App and command line applying different rules to the same action.
- Any regression in the app for users who never touch the command line.
- Auto-accepting something a careful curator would have stopped on.

## Constraints

- No new package dependencies. `cli`, `processx`, and `jsonlite` are already
  imported.
- Supported Nextflow range and version pinning unchanged; the launcher must
  go through the existing checks.
- Nextflow must start from the project directory (relative database path).
- Two writers on one project database must be avoided.
- No non-ASCII characters in code.

## Approaches considered

1. **Shared core (chosen).** Extract the actions the command line needs from
   app server code into plain functions; app observers and commands both call
   them. Deep editing stays in the app. One set of rules; costs careful app
   regression testing.
2. **Parallel terminal code, app untouched.** Zero app risk, but duplicates
   rules such as lock eligibility, which then drift. Rejected.
3. **Full separation of all app actions.** Cleanest end state, far too large,
   and covers edits that belong in the app anyway. Rejected.

## How users get the command

The command is a small shell script that runs
`Rscript -e 'MitoPilot:::cli_main()' "$@"` with a specific Rscript.

- **Container-based install** (R and the MitoPilot package on the host,
  containers used per pipeline step; today's standard install). The script
  ships inside the package. One-time setup after installing the package:
  `Rscript -e 'MitoPilot::install_cli()'`. This writes `mitopilot` into
  `~/.local/bin` (or `dir =`), with the path of the Rscript that ran it
  baked in, and warns if that directory is not on PATH. Re-run after moving
  R. Package updates need no re-run, since the script only hands off to the
  installed package.
- **Native (no-container) install.** `install_mitopilot_native.sh` writes
  `mitopilot` into `<prefix>/envs/mitopilot/bin`. After
  `source <prefix>/activate.sh` the command is on PATH. No extra step.
- **Docker image.** `mitopilot` is also placed on PATH in the image so
  `singularity exec <image> mitopilot ...` works for the site setups that
  already run R that way. Not a supported way to launch pipeline runs.

In all cases the command acts on the project in the current directory, or
the one given by `--project DIR`.

## Command surface

```
mitopilot init [DIR] [--mapping FILE --data DIR | --assemblies ... | --annotations ...]
               [--test reads|assemblies|annotations] [--n N] [--full-size] [--no-raw-data]
               [--executor X] [--container X] [--force] [--no-metadata]
               [--settings settings.json]
mitopilot status [--json] [--stage assemble|annotate|export] [--parked]
mitopilot run assemble|annotate [--no-resume] [--submit] [--dry-run]
mitopilot accept SAMPLE... [--stage ...] [--override [--reason TEXT]]
mitopilot accept --all-parked --override [--stage ...]
mitopilot hold|requeue SAMPLE... [--stage ...]
mitopilot rules [--min-coverage N] [--accept-linear | --no-accept-linear]
mitopilot group NAME --samples ... | --all-accepted
mitopilot metadata fetch [--source geome|gbif|ncbi]
mitopilot export GROUP [--flagged ask|stop|drop|keep]
mitopilot auto [--group NAME] [--submit] [--dry-run]   # run, auto-accept, run, auto-accept, export
mitopilot cluster setup | list | show NAME | edit NAME | check [NAME] | apply NAME
mitopilot opts list | show STAGE SET | set STAGE SET --key value ... |
               assign STAGE SET --samples ... | edit [STAGE] |
               export [FILE] | apply FILE
mitopilot review [--stage ...] [--samples ...]
mitopilot app [--host H] [--port P] [--no-browser]
```

- `--test` maps to `new_test_project`, `new_test_project_userAsmb`,
  `new_test_project_userAnnot`; other flags pass through.
- Default output is a readable table; `--json` on `status` and `opts show`.
- Exit codes: 0 done and nothing parked; 1 error; 2 finished but samples
  are parked. Scripts branch on 2.
- Any change that makes finished samples rerun prints the count and needs
  confirmation, or `--yes` in scripts.
- Flag reading is home-made (`--key value`, `--flag`, `--no-flag`). `opts
  set` flags are generated from the settings table columns.

## Pass rules (auto-accept)

Evaluated by one function per stage that returns, per unit, pass or a list
of reasons. `status` shows the reasons; `auto` and `run --auto-accept` act
on them.

- **Assembly:** run succeeded; one path; one scaffold; circular; mean
  coverage of the assembly at or above a per-project threshold (new setting,
  default 10x, matching the existing path-scoring warning); reference match
  not flagged poor.
- **Annotation:** validation passed and zero curation warnings. Any warning
  of any kind parks the sample. `status` also prints a tally of warning
  types across the project, for future tuning.
- **Export group:** no outlier flags; no internal stops (see Export).

Fixed conservative defaults. Per-project on/off switches for individual
rules, stored in the project; first switch: accept a clean linear assembly
and mark it partial (default off, so linear assemblies are parked).

A pass results in the same lock the app applies. Auto-accepted units are
recorded as such so `status` and the app can tell them from human accepts.

### Overrides

A parked or flagged unit never moves forward on its own. It moves only by an
explicit human override, either:

- from the command line: `accept SAMPLE --override`, or in bulk with
  `accept --all-parked --override`; or
- per sample in the `review` walk-through ("accept anyway"), or by fixing
  and locking it in the app.

Overrides are recorded (time, optional reason, which rules were failing) and
shown in `status`.

## Export

Export choices that the app's export window holds only in memory today
(start and stop offsets, identity percent, alignments on or off, which
genes) become a saved export settings stage, read by both the app window
and the command line, and covered by `opts`.

`export GROUP` runs the outlier check first and always writes the HTML
alignment report, printing its path. If nothing is flagged it exports.
If samples are flagged:

- At a terminal (`--flagged ask`, the default there) it lists the flagged
  samples and genes, points to the HTML report, and asks:
  1. review in app: opens the app with the outlier review for that group
     open; when the window closes the check reruns and the prompt repeats
     if anything is still flagged
  2. drop flagged samples and export the rest
  3. keep all and export (an override, recorded as above)
  4. cancel
- In a script (no terminal) the default is `--flagged stop`: nothing is
  exported, exit code 2. `--flagged drop` and `--flagged keep` are the
  explicit scripted choices. `auto` follows the same rule.

`metadata fetch` wraps the existing GEOME, GBIF, and NCBI fetchers. Mapping
source fields to columns stays in the app.

## Clusters

`run` and `auto` always start Nextflow directly in the current process,
whatever the machine. They never poll a scheduler. On a cluster the user
runs them inside a scheduler job; `--submit` writes that job script and
submits it, and `--dry-run` prints the script and the Nextflow command
without submitting. Nextflow then hands pipeline steps to the scheduler
through the project `.config`, as today. `auto --submit` submits one job
that runs the whole chain, so it survives logout. Without `--submit` on a
machine where a scheduler is detected, `run` warns that it is about to run
the pipeline head process on the current node and continues.

Two things must be right per cluster: the Nextflow executor settings
(profile, already handled by `generate_config()`, saved profiles, and the
bundled slurm, sge, pbs, lsf, Hydra, and SEDNA templates) and the header of
the submit script (already handled per project by
`.mitopilot_submit.template`). `mitopilot cluster` puts both behind one
command and stores the submit template with the saved profile, so `init
--executor NAME` gives a project both.

- `cluster setup`: guided questions (scheduler type, queue, account or
  project code, extra scheduler options, parallel environment for SGE,
  container engine and cache folder or native install folder; then head-job
  queue, time limit, memory, and any module or `source` lines needed before
  the run). Saves a named profile. Ends by offering `cluster check`.
- `cluster list`, `cluster show NAME`.
- `cluster edit NAME`: opens the saved Nextflow config or submit template in
  `$EDITOR`. This is the escape hatch for clusters the questions do not
  cover; both files are plain text and fully user-owned.
- `cluster check [NAME]`: staged smoke test that stops at the first failure
  and prints the scheduler's own error text with a hint:
  1. scheduler commands, Nextflow (supported version), Java, and container
     engine or native environment found
  2. config parses (`nextflow config`)
  3. submit script accepted: submits a job that only prints the hostname
  4. Nextflow can hand a step to the scheduler: a one-task trivial pipeline
     using the profile, run in the same environment a real step gets
     (container or native), calling `Rscript -e 'library(MitoPilot)'`
  Takes minutes, so a bad queue name or missing bind path is found before a
  real run.
- `cluster apply NAME`: switch an existing project to a saved profile
  (existing config migration).

The app's "Write Submission Script" path keeps working and reads the same
templates.

## Shared core (extracted from app code)

Plain functions taking a project path or connection, returning the rows they
wrote. The app observers are changed to call them and then refresh their
reactive copies.

- Launch a run: queue count, executor and native environment, version pin,
  start Nextflow or build and submit the scheduler script, write the run
  report. Progress in the terminal uses the existing log-based progress
  functions.
- Set sample state; lock and unlock for assembly and annotation, including
  existing blocks (missing output directory, more than one path locked).
- Save a settings set; assign a set to samples (including requeue).
- Assign export group; pre-export checks.
- Pass rule evaluation.

Everything else in the app is left as is.

## Settings

Three routes onto the same save functions:

- **Guided editor** `opts edit`: numbered menus, stage then set then
  setting. Shows current value, default, allowed choices; `?` shows the
  explanation, taken from the existing in-app help text. Settings that only
  apply to a given choice are hidden until relevant. Nothing is saved until
  confirmed. On save it prints the equivalent one-line command. Refuses to
  run without a terminal and points to the one-line commands.
- **One-line commands** `opts set / assign / list / show`.
- **File** `opts export / apply`: JSON holding every set and the sample
  assignments; `apply` shows the changes first. `init --settings` starts a
  project from it. Sample assignment is also accepted as extra columns in
  the mapping file at `init`.

Curation: choosing a rule set fills its default thresholds; individual
thresholds can be overridden. Invalid names or values stop with the valid
choices listed and change nothing.

## Review walk-through

`mitopilot review` steps through parked units in pipeline order. For each it
prints sample, stage, and reasons, then a menu:

- open in app
- terminal-only choices where sensible: requeue with a different settings
  set, hold, accept anyway, skip
- quit

"Open in app" starts the app in the foreground on the matching tab, table
filtered to that sample (filter clearable), details window open. The app
already stops when the last browser window closes; the command then reports
what changed on that sample and moves on. It refuses to open the app while
a pipeline run is active on the project.

### Local and remote

- Desktop session: open the browser.
- SSH session (detected from the environment): bind for tunnelling, print
  the tunnel command and link using the existing helper, wait.
- `--host`, `--port`, `--no-browser` override. The walk-through reuses one
  port for its whole run so the tunnel and link stay valid across samples.
- Site setups (Hydra, SEDNA) go through the existing handling.

### App behaviour when started from the command line

| Feature | Focused review launch | `mitopilot app` |
|---|---|---|
| Update (run workflows) | disabled, with a note to use `mitopilot run` | disabled |
| Accept, lock, path pick, consensus, join, gene edits, flags, notes | open | open |
| Settings sets, set sample to ready | open; command reports requeued count afterwards | open |
| Liftover | open | open |
| Export button | disabled | open |
| Outlier review, run reports, file browser, metadata | open | open |

`MitoPilot()` started the normal way is unchanged.

## Testing

- Unit tests for every extracted function, pass rules, flag reader, and
  settings validation.
- End-to-end command-line run on the three test projects in a scratch copy:
  `init --test`, `auto`, exit codes, parked reasons.
- **App regression, required for every extracted action:** open the real app
  in a browser on a scratch project copy and click through the action before
  and after (state change, lock and unlock with their blocks, each settings
  window save and assign, export group, run launch, export), comparing the
  database rows written. Both the read-based and user-assembly apps.
- Focused launch, disabled buttons, close-window hand-back, and the remote
  path (tunnel on a second machine or loopback simulation).
- `cluster check` and `--submit` on at least one real scheduler (Hydra is
  available); other scheduler types by `--dry-run` script comparison against
  the app's existing output.

## Build order

1. Command, `init`, `status`, `run`, pass rules, `accept` and friends,
   `group`, `metadata fetch`, `export` (scripted `--flagged` choices and
   the terminal prompt without the app option), `auto`, `cluster`.
2. `opts` (commands, file, guided editor).
3. `review`, focused app launch, `mitopilot app`, button restrictions,
   local and remote handling, and the "review in app" choice in `export`.

Each stage is usable on its own. Between 1 and 3, parked samples are handled
by opening the app the normal way.

## Left to the builder

- Internal function and file names; table or column used to record
  auto-accept and rule switches (with a `backwards_compatibility` migration).
- Exact wording of menus, tables, and messages.
- How the focused sample and command-line mode are passed into the app
  (existing launch options are available).
- How a pipeline run in progress is detected.
- Detailed rule checks per clade where curation warnings differ.
