# POWHEG

ALICE copy of POWHEG. Please refer to http://powhegbox.mib.infn.it/ for the official sources.

# importing the sources

The upstream sources are maintained in an svn repository and are
divided into a core part and user processes. In order to update
POWHEG, replace the content of the main directory with the one
obtained by:

```
svn export --force --username anonymous --password anonymous svn://powhegbox.mib.infn.it/trunk/POWHEG-BOX-V2 .
svn export --force --username anonymous --password anonymous svn://powhegbox.mib.infn.it/trunk/User-Processes-V2/hvq
svn export --force --username anonymous --password anonymous svn://powhegbox.mib.infn.it/trunk/User-Processes-V2/W
svn export --force --username anonymous --password anonymous svn://powhegbox.mib.infn.it/trunk/User-Processes-V2/Z
svn export --force --username anonymous --password anonymous svn://powhegbox.mib.infn.it/trunk/User-Processes-V2/dijet
```
and tag the version with the svn revision.

In case you need to add a patch on top of this, create a branch with:

```
git checkout -b alice/<version> <version>
```
