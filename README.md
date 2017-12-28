## Happah Core

### Contributing

New developers who wish to contribute to the project can get started by executing the following commands on an Ubuntu 16.04.2 machine:

```
sudo apt install git dh-autoreconf libglm-dev libeigen3-dev liblpsolve55-dev libboost-all-dev autoconf-archive
mkdir -p ${HOME}/Workspace/include
cd ${HOME}/Workspace
wget -O boost_1_64_0.tar.gz https://dl.bintray.com/boostorg/release/1.64.0/source/boost_1_64_0.tar.gz
tar -C include --strip=1 -xzvf boost_1_64_0.tar.gz boost_1_64_0/boost/spirit boost_1_64_0/boost/spirit.hpp
git clone http://github.com/happah-graphics/happah-core.git
cd happah-core
./bootstrap
./configure --prefix=${HOME}/Workspace
export CPATH="${HOME}/Workspace/include"
export LD_LIBRARY_PATH="${HOME}/Workspace/lib"
make install
```

Make your changes and run ``` make ``` to compile the library and ``` make install ``` to install the library into the include and lib directories.  Finally, run ``` git push origin master ``` to upload your changes to Github.

If you have a release-ready version, tag it by executing ``` git tag -a v0.1 -m "version 0.1" ``` and upload the tag to Github using ``` git push origin v0.1 ``` to push a specific tag or ``` git push origin --tags ``` to push all tags at once.

