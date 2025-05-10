# uyduTasarimi
This repository hosts a thermal simulation code for satellite systems

# Use

We need `helezon` to run the thermal simulations in this project.

## Conda

This library can be installed directly with conda. We can make a conda environment for fenicsx and then install this library into it.

```bash
conda create -n uydu python=3.13.0
conda activate uydu
conda install pip
conda install -c conda-forge fenics-dolfinx=0.9.0 pyvista=0.44.1 # Linux and macOS
# conda install -c conda-forge fenics-dolfinx pyvista pyamg # Windows
```

We need `helezon` to run the thermal simulations in this project. Go to helezon folder and run

```bash
pip3 install -e .
```

If you receive the error

```bash
ImportError: /usr/lib/x86_64-linux-gnu/libstdc++.so.6: version `GLIBCXX_3.4.31' not found
```

You need to to 
```bash
strings $CONDA_PREFIX/lib/libstdc++.so.6 | grep GLIBCXX_3.4.31
```

then 

```bash
mkdir -p $CONDA_PREFIX/etc/conda/activate.d
echo "export LD_PRELOAD=\$CONDA_PREFIX/lib/libstdc++.so.6" > $CONDA_PREFIX/etc/conda/activate.d/libstdcpp_preload.sh

mkdir -p $CONDA_PREFIX/etc/conda/deactivate.d
echo "unset LD_PRELOAD" > $CONDA_PREFIX/etc/conda/deactivate.d/libstdcpp_preload.sh

```
