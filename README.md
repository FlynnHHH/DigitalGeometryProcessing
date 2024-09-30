# Tutte Parameterization

## Dependency

* [Eigen](http://eigen.tuxfamily.org/)

## Usage

```SHELL
git clone https://github.com/FlynnHHH/TutteParameterization.git
```

Edit lines **9** of CmakeLists.txt to set the values of **EIGEN_PATH**

```SHELL
mkdir build && cd build
cmake ../
make
./TUTTE_PARAMETERIZATION path/to/input_mesh.obj /path/to/output_mesh.obj
```
