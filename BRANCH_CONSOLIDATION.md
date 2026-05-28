# Consolidating Multiple Branches into One Main Pull Request

This guide explains how to move work from several feature branches into one clean branch so you can open **one main pull request**.

> In this repository snapshot, the current local branch is `work`, and it already contains the converter, README, sample inputs, and example outputs. If your GitHub account shows additional branches, follow the steps below on your Windows computer where those branches exist.

## Goal

Create one branch, for example `main-iscn-converter-pr`, that contains all changes you want reviewed. Then open one pull request from that branch.

## Before you start

Open PowerShell in your repository folder:

```powershell
cd "C:\Users\YOUR_NAME\Documents\ISCN-Conversion-"
```

Check where you are:

```powershell
git status
git branch --all
```

If `git status` shows uncommitted changes, either commit them or save them before switching branches:

```powershell
git add .
git commit -m "Save current work before branch consolidation"
```

## Option A: Easiest approach if one branch already has everything

If one branch already contains all the files you want in the pull request, use that branch as the PR branch.

```powershell
git switch branch-that-has-everything
git pull
```

Then push it:

```powershell
git push -u origin branch-that-has-everything
```

Open one pull request from `branch-that-has-everything` into `main`.

## Option B: Combine several branches with cherry-pick

Use this when each branch has one or more useful commits that you want to combine.

### Step 1: Update your local repository

```powershell
git fetch --all --prune
```

### Step 2: Create one new PR branch from `main`

```powershell
git switch main
git pull origin main
git switch -c main-iscn-converter-pr
```

### Step 3: Find the useful commits on each old branch

```powershell
git log --oneline main..old-branch-name
```

Repeat that command for each branch you want to merge.

### Step 4: Cherry-pick the useful commits

Use the commit hashes from the previous step:

```powershell
git cherry-pick COMMIT_HASH_1
git cherry-pick COMMIT_HASH_2
git cherry-pick COMMIT_HASH_3
```

If there is a conflict:

1. Open the conflicted files.
2. Keep the version you want.
3. Remove conflict markers like `<<<<<<<`, `=======`, and `>>>>>>>`.
4. Run:

```powershell
git add .
git cherry-pick --continue
```

If you picked the wrong commit and want to stop:

```powershell
git cherry-pick --abort
```

### Step 5: Run the converter checks

```powershell
py .\iscn_to_plink_cnv.py --run-tests
py .\iscn_to_plink_cnv.py --raw .\sample_raw.tsv --cytoband .\sample_cytoband.tsv --outdir .\out
```

### Step 6: Push the one consolidated branch

```powershell
git push -u origin main-iscn-converter-pr
```

Then open one pull request from `main-iscn-converter-pr` into `main`.

## Option C: Combine whole branches with merge

Use this if you want everything from each branch, not selected commits.

```powershell
git fetch --all --prune
git switch main
git pull origin main
git switch -c main-iscn-converter-pr
git merge old-branch-1
git merge old-branch-2
git merge old-branch-3
git push -u origin main-iscn-converter-pr
```

If conflicts happen, resolve them, then run:

```powershell
git add .
git commit
```

## Closing extra pull requests

After the one consolidated PR is open:

1. Go to GitHub.
2. Open each older pull request.
3. Leave a comment such as: `Closing this PR because the changes were consolidated into #NEW_PR_NUMBER.`
4. Close the older pull request.

Do not delete old branches until the new consolidated pull request has been reviewed and merged.

## Recommended path for this project

For this converter project, the safest path is usually:

```powershell
git fetch --all --prune
git switch main
git pull origin main
git switch -c main-iscn-converter-pr
```

Then cherry-pick only the commits that add or improve:

- `iscn_to_plink_cnv.py`
- `README.md`
- `sample_raw.tsv`
- `sample_cytoband.tsv`
- `out/plink_cnv_output.tsv`
- `out/review_flags.tsv`
- `out/optional_debug_extracted_events.tsv`

This creates one clean PR and avoids accidentally bringing in unrelated work.
