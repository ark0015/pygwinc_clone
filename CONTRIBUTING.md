# Contributions to pygwinc

The pygwinc project welcomes your contributions.  Our policy is that
all contributions should be peer-reviewed.  To facilitate the review,
please do not commit your work to this repository directly.  Instead,
please fork the repository and create a [merge
request](https://git.ligo.org/gwinc/pygwinc/-/merge_requests/new)
against the main pygwinc master branch.

When submitting code for merge, please follow good coding practice.
Respect the existing coding style, which for `pygwinc` is standard
[PEP8](https://www.python.org/dev/peps/pep-0008/) (with some
exceptions).  Make individual commits as logically distinct and atomic
as possible, and provide a complete, descriptive log of the changes
(including a top summary line).  Review efficiency follows code
legibility.

`pygwinc` comes with a validation command that can compare budgets
between different versions of the code and produce a graphical report.
`pygwinc` stores a reference to the git commit that is considered the
latest good reference for the budget calculation functions, and the
validation is run against that reference by default.  The command can
be run with:
```shell
$ python3 -m gwinc.test
```
Use the '--plot' or '--report' options to produce visual comparisons
of the noise differences.  The comparison can be done against an
arbitrary commit using the '--git-ref' option.  Traces for referenced
commits are cached, which speeds up subsequent comparison runs
significantly.

Commits for code changes that modify budget calculations should be
followed up updates to test reference. This can be done with the
'--update-ref' option to the test command:
```shell
$ python3 -m gwinc.test --update-ref
``` 
If no specific reference is provided, a pointer to the last commit
will be made.  *Note: reference updates need to be made in a separate
commit after the commits that change the noise calculations* (git
commits can't have foreknowledge of their own commit hash).

Once you submit your merge request, an automated pipeline will run the
test command to validate your changes.  If there are budgets
differences but the reference has not been updated, the pipeline will
fail and indicate that a reference update is required.  Updates to the
reference will also cause a validation pipeline failure, but these
failures can be resolved through reviewer approval of the changes (see
below).


## For admins: approving noise curve changes

As discussed above, merge requests that generate noise changes will
cause pipeline failures.  If the MR properly includes a reference
update, then the failure should only be in the approval check.  A
report will be generated comparing all noise changes against the
target branch (usually 'master').  See the 'View exposed artifacts'
menu item in the pipeline report.  Once you have reviewed the report
and the code, and understand and accept the noise changes, click the
'Approve' button in the MR.  Once sufficient approval has been given,
re-run the `review:check_approval` pipeline job, which should now pick
up that approval has been given and allow the pipeline to succeed.
Once the pipeline succeeds the merge can be enacted.  Click the
'Merge' button to finally merge the code.
