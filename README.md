# CMake Wrapper of implementation "Isosurfaces Over Simplicial Partitions of Multiresolution Grids"

This is a cmake wrapper of Josiah Manson's source code of the paper "Isosurfaces Over Simplicial Partitions of Multiresolution Grids". This method is super nice and solve almost every desired properties of iso-surface generation: 

- Intersection-free
- Manifold
- Sharp Feature preservation.

## Dependencies
I use [libigl](https://libigl.github.io/) for mesh processing, [polyscope](https://github.com/nmwsharp/polyscope.git) for visualization and [CLI11](https://github.com/CLIUtils/CLI11) for command line arguments parsing. All the dependencies are built through FetchContent. No extra action is needed. 

## Compile
Assume you are use linux-like machine, compile this project using the standard cmake routine:

    git clone https://github.com/zhenchen-jay/isosurface-wrapper.git
    mkdir build
    cd build
    cmake ..
    make -j4
For windows, you can use CMake Gui + Visual Studio to achieve this.

## Run
```
./bin/isosurface_extract_bin --csg
```
After that, you will see a polyscoped rendered mesh corresponding to the figure 3(b) in the paper.
If use wish to run more general model, like the isosurface of any input mesh (.obj), I also implemented a fake CSG subclass call GeneralShape. You can try with 
```
./bin/isosurface_extract_bin -i [input obj path] -d [offset distance, w.r.t. diagnoal size of bounding box] --minDepth (optional, by default, it is 4) --maxDepth (optional, by default, it is 7) -o [output file, optional, by default it is input_path_minDepth_maxDepth_offset.obj]
```

## Acknowledgement
Thanks Josiah Manson for providing his implementation of the paper. You can check the original implementation here: http://josiahmanson.com/research/iso_simplicial/.
