#include <Halide.h>
#include <mpi.h>
int main(int argc, char **argv) {
// We’ll take the pipeline from the previous example and
// turn it into a distributed pipeline using the new
// language features.
// This is now an MPI program, so we need to initialize
// the MPI library.
MPI_Init(&argc, &argv);
// Declare and allocate a 100x100 input and output
// image. Now that these are DistributedImages, the
// 100x100 refers to the global extents, not the extents
// on each rank.
Halide::DistributedImage<int> input(100, 100),
output(100, 100);
// The pipeline definition (i.e. the algorithm) is
// identical to the non-distributed example. Nothing
// changes here.
const float factor = 1.25f;
Halide::Func brighten;
Halide::Var x, y;
brighten(x, y) = factor * input(x, y);
// Define a second pipeline stage to convert the new
// pixel values back to integers and clamp their values
// to [0,255].
Halide::Func clampint;
clampint(x, y) = clamp(Halide::cast<int>(brighten(x, y)),
0, 255);
// Now distribute the pipeline on the y dimension. Each
// rank will be responsible for computing a contiguous
// "slice" of rows for the ’clampint’ stage. The
// ’brighten’ stage is unscheduled, which means the
// ’brighten’ function will be inlined into ’clampint’.
clampint.distribute(y);
// The definition and scheduling of the distributed
// Halide pipeline is complete. Now we can specify the
// data distribution. This has to occur after scheduling
// the pipeline so that the input image can be allocated
// with room for any border exchanges (here there are
// none). Here we specify both the input and output
// images should also be distributed on the y dimension.
input.set_domain(x, y);
input.placement().distribute(y);
input.allocate();
output.set_domain(x, y);
output.placement().distribute(y);
output.allocate();
// Compile.
clampint.compile_jit();
// Initialize the pixel values with an arithmetic
// progression. Now the x and y iterate over the
// rank-local section of ’input’. We can access the
// global extents with ’input.global_height()’ and
// ’input.global_width()’.
for (int y = 0; y < input.height(); y++) {
for (int x = 0; x < input.width(); x++) {
// Because x and y are local coordinates, we use
// the DistributedImage::global() function to
// convert local x to global x and local y to
// global y.
const int global_x = input.global(0, x);
const int global_y = input.global(1, y);
input(x, y) = global_x + global_y;
}
}
// Execute it into the rank-local portion of the output
// image.
clampint.realize(output);
MPI_Finalize();
return 0;
}
