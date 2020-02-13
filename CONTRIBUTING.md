# Contributions to pygwinc

The pygwinc project welcomes your contributions.  Our policy is that
all contributions should be peer-reviewed.  To facilitate the review,
please do not commit your work to this repository yourself.  Instead,
fork the repository and send a merge request.

`pygwinc` comes with a test suite that will calculate all canonical
IFO budgets with the current state of the code, and compare them
against cached hdf5 traces (stored in the repo with
[git-lfs](https://git-lfs.github.com/) ( [see
also](https://wiki.ligo.org/Computing/GitLFS)).  Contributors should
run the tests before submitting any merge requests:
```shell
$ python3 -m gwinc.test
```
The test suite will be run automatically upon push to git.ligo.org as
part of the continuous integration, and code that fails the CI will
not be merged, unless failures are well justified.

If your change affects the shape of a noise curve, your commit message
should make note of that, and provide a justification.  It will also
be necessary to update the reference curves stored in cache, but that
should be done in a separate commit not a part of the original merge
request, so that all reviewers can see the changes introduced in the
CI test failure report.
