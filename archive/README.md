# archive/

The 43 original R scripts, kept byte-for-byte as they were written.

**They are not part of the pipeline.** They are kept so that every result in the
paper can be traced back to the exact code that produced it, including the
exploratory detours. The maintained, re-runnable versions are the numbered
modules in `R/`.

Nothing here is executed by `run_all.R`, and nothing here is needed to
reproduce a figure. Several scripts do not run as-is (hard-coded paths on
`D:/`, `E:/`, `F:/`, `G:/`, undefined objects, one syntax error) — that is
expected, and it is why the pipeline was rebuilt.

`docs/migration_notes.md` maps each of these files to the module that replaced
it and lists the defects that were fixed.
