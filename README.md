## Happah Core

### Contributing

New developers who wish to contribute to the project can get started by executing the following commands on an Ubuntu machine:

```
sudo apt-get install git devscripts equivs
git clone http://github.com/happah-graphics/happah-core.git
cd happah-core
git checkout -b package/wily origin/package/wily
sudo mk-build-deps -i
git checkout master
./bootstrap
./configure or ./configure --prefix=${dir}
make && make install
```

Make your changes and run ``` make && make install ``` to update the library.  Finally, run ``` git push origin master ``` to upload your changes to Github.

If you have a release-ready version, tag it by executing ``` git tag -a v0.1 -m "version 0.1" ``` and upload the tag to Github using ``` git push origin v0.1 ``` to push a specific tag or ``` git push origin --tags ``` to push all tags at once.

To start programming, save the following code in a text file called main.cpp

```
#include <iostream>
#include <happah/io/writers/WriterOFF.h>
#include <happah/geometries/TriangleMesh.h>

int main() {
     std::cout << "INFO: Printing triangle mesh.\n";
     happah::TriangleMesh3D mesh({{{0, 0, 0}}, {{1, 0, 0}}, {{0, 1, 0}}}, {0, 1, 2});
     happah::WriterOFF::write(mesh, "single-triangle.off");
     std::cout << "INFO: Done printing triangle mesh.\n";
     return 0;
}
```

and compile it by executing ``` g++ -o a main.cpp -I/usr/local/include -fcilkplus -std=c++1y ```.

### Building the Debian Package

Package maintainers can build the Debian package by executing the following commands on an Ubuntu machine:

```
sudo apt-get install git git-buildpackage devscripts equivs
git clone http://github.com/happah-graphics/happah-core.git
cd happah-core
git checkout -b package/wily origin/package/wily
sudo mk-build-deps -i
gbp buildpackage -us -uc
```

To clean the directory of all build files, run ``` dh clean ```.

To update the package to a new version of the library, execute:

```
git checkout package/wily
git merge v0.1 --squash
gbp buildpackage -us -uc
```

Once the package builds successfully, run ``` gbp buildpackage -us -uc --git-tag ``` to tag the package version.  ``` --squash ``` prevents git from copying the commit log history into the package/wily branch.

To create the build-deps package, execute:

```
sudo apt-get install devscripts equivs
mk-build-deps happah-core_0.1-1.dsc
```

Then, upload the resulting package happah-core-build-deps_0.1-1_all.deb into the dropbox branch.

Alternatively, instead of building and distributing the build-deps package, simply run ``` sudo mk-build-deps -i ``` while in the package/wily branch.

