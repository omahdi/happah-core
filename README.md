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

To start programming, save the following code in a text file called main.cpp

```
#include <iostream>
#include <happah/format.h>
#include <happah/geometries/TriangleMesh.h>

int main() {
     std::cout << "INFO: Writing triangle mesh.\n";
     auto mesh0 = happah::make_triangle_mesh<VertexP3>({{{0, 0, 0}}, {{1, 0, 0}}, {{0, 1, 0}}}, {0, 1, 2});
     happah::format::off::write(mesh0, "single-triangle.off");
     std::cout << "INFO: Done writing triangle mesh.\n";
     std::cout << "INFO: Reading triangle mesh.\n";
     auto content = happah::format::off::read("single-triangle.off");
     auto mesh1 = happah::make_triangle_mesh<VertexP3>(content);
     std::cout << mesh1 << '\n';
     std::cout << "INFO: Done reading triangle mesh.\n";
     return 0;
}
```

and compile it by executing ``` g++ main.cpp -std=c++1y -I${HOME}/Workspace/include -L${HOME}/Workspace/lib -lhappah -lboost_iostreams ```.  Then, execute ``` ./a.out ```.

