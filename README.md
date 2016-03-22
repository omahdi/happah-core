## Happah Core

### Contributing

New developers who wish to contribute to the project can get started by executing the following commands on an Ubuntu machine:

```
sudo apt-get install git
wget https://github.com/happah-graphics/happah-core/raw/dropbox/happah-core-build-deps_0.1-1_all.deb
sudo dpkg -i happah-core-build-deps_0.1-1_all.deb
sudo apt-get -f install
git clone http://github.com/happah-graphics/happah-core.git
cd happah-core
./bootstrap
./configure or ./configure --prefix=${dir}
make && make install
```

Make your changes and run ``` make && make install ``` to update the library.  Finally, run ``` git push origin master ``` to upload your changes to Github.

If you have release ready version, tag it by executing ``` git tag -a v0.1 -m "version 0.1" ``` and uploading the tag to Github using ``` git push origin v0.1 ``` to push a specific tag or ``` git push origin --tags ``` to push all tags at once.

### Building the Debian Package

Package maintainers can build the Debian package by executing the following commands on an Ubuntu machine:

```
sudo apt-get install git git-buildpackage
wget https://github.com/happah-graphics/happah-core/raw/dropbox/happah-core-build-deps_0.1-1_all.deb
sudo dpkg -i happah-core-build-deps_0.1-1_all.deb
sudo apt-get -f install
git clone http://github.com/happah-graphics/happah-core.git
cd happah-core
git checkout -b package/wily origin/package/wily
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

