# Replace the GitHub repository safely

Use this when the GitHub web editor or manual upload has corrupted line breaks.

## One-command method

From the unzipped workflow folder:

```powershell
powershell -ExecutionPolicy Bypass -File .\scripts\replace_github_repo.ps1
```

The script will:

1. clone `https://github.com/benchtopc/Dual-Barcode-Demux.git`
2. preserve `.git`
3. delete the old broken files
4. copy this clean workflow into the repo root
5. commit
6. push to `main`

## Dry run

```powershell
powershell -ExecutionPolicy Bypass -File .\scripts\replace_github_repo.ps1 -DryRun
```

## After upload

Import into EPI2ME with:

```text
https://github.com/benchtopc/Dual-Barcode-Demux
```

Do not use the `.git` suffix in EPI2ME.
