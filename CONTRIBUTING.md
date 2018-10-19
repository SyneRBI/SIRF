Contributing
============

Please help us by finding problems, discussing on the mailing lists, contributing documentation,
bug fixes or even features. Below are some brief guidelines.

## Reporting a problem

Please use our [issue-tracker]: https://github.com/CCPPETMR/SIRF/issues

## Submitting a patch

For contributing any code or documentation that is non-trivial, we require a
signed Contributor License Agreement, stating clearly that your
conributions are licensed appropriately. This will normally need to be signed by your
employer/university, unless you own your own copyright.
You will have to do this only once. Please contact us for more information.

Please keep a patch focused on a single issue/feature. This is important to keep our history clean,
but will also help reviewing things and therefore speed-up acceptance.

### Process

This is our recommended process. If it sounds too daunting, ask for help.

1. Create a new issue (see above). State that you will contribute a fix if you intend to do so.
2. Create a [fork](https://help.github.com/articles/fork-a-repo) on github and work from there.
3. Create a branch in your fork with a descriptive name and put your fixes there. If your fix is
simple you could do it on github by editing a file, otherwise clone your project (or add a remote
to your current git clone) and work as usual.
4. Use [well-formed commit messages](http://tbaggery.com/2008/04/19/a-note-about-git-commit-messages.html)
for each change (in particular with a single "subject" line
followed by an empty line and then more details).
5. Push the commits to your fork and submit a [pull request (PR)](https://help.github.com/articles/creating-a-pull-request)
(enable changes by project admins.) Be prepared to add further commits to your branch after discussion.
In the description of the PR, add a statement about which Issue this applies to
using [a phrase such that github auto-closes the issue when merged to master](https://help.github.com/articles/closing-issues-using-keywords/).
6. After acceptance of your PR, go home with a nice warm feeling.

## Project rules

- Only one official, stable, up-to-date branch: **master**
    + Essentially "latest stable beta version with no known bugs
      since the last official release version"
    + Never knowingly add a bug to **master**
- Any work-in-progress commits should be in their own branches.
- GitHub assigns a unique number to each issue, c.f. the [issue-tracker].
- A pull request (PR) is an issue with an associated branch,
  c.f. [pull-requests]. Even for "internal" development, we prefer a PR for
  a branch to allow review and discussion.
- Branches and PRs are kept small (ideally one 'feature' only) and branch from **master**,
  not from another branch, unless required. This allows
  commenting/improving/merging this branch/PR
  independent of other developments.
- Discussions on issues and PRs are forwarded to the
  <CCP-PETMR-DEVEL@jiscmail.ac.uk> mailing list daily.
    + Forwarded from github via the [googlegroup],
      which is also a backup in case github dies.
- Contributions of new features should also update documentation and release notes. After version 1.0,
  this needs documentation needs to state something like "introduced after version 1.xxx".
- We prefer issues to be opened via [github][issue-tracker] due to the following reasons:
    + Ensures issues will never get lost in emails
        * Facilitates issue status tracking
    + Allows focused comments/discussion
        * Easy cross-referencing of related issues, PRs, and commits
    + The mailing list gets notified within 24 hours.

[issue-tracker]: https://github.com/CCPPETMR/SIRF/issues
[pull-requests]: https://github.com/CCPPETMR/SIRF/pulls
[googlegroup]: https://groups.google.com/forum/#!forum/ccp-petmr-codebot


## Submodules

> [...] A submodule is a repository embedded inside another repository. The
> submodule has its own history [...]

This project uses [submodules] to point to infrequently changed/large datasets.

1. How to pull updates to SIRF
```bash
# first time
SIRF$ git clone https://github.com/CCPPETMR/SIRF --recursive
# subsequently
SIRF$ git pull
SIRF$ git submodule update --init --recursive
# or
SIRF$ git pull --recurse-submodules
```
2. How to push updates to the data submodule in SIRF
```bash
SIRF$ cd data
SIRF/data$ # create/change some files
SIRF/data$ git commit -m "add some data" && git push && cd ..
SIRF$ git add data && git commit -m "update submodule" && git push
```
3. What happens when switching branches
    + When switching branches in SIRF, the data submodule will not be touched.
      `git submodule update` needs to be manually run in order to ensure the
      correct submodule commit is also checked out
	+ WARNING: This is why `git add --all` is bad practice (if you forget to
	  update submodules when needed, then committing all modifications will
	  include the current wrong submodule commit)

[submodules]: https://git-scm.com/docs/gitsubmodules
