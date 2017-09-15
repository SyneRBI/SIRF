Contributing
============

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
