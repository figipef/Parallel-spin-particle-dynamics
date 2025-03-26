# -Parallel-spin-particle-dynamics


# Quick Input laser tutorial

Type 0 - Constant fields (intialize both direct E and B or through E and k vector)

Type 1 - Fields with a certain frequency (intialize both direct E and B and frequency or through E and k vector)

Type 2 - Continuous Laser (initialize through E and k vector)

Type 3 - Laser in a packet (initialize through k vector and a pakcet frequency and length)

# Particle Creation tutorial

-> Spin

Type 0 - Uniform random direction with unit norm

Type 1 - Specific direction chosen from input (normalized inside the code)

Type 2 - Von Miser-Fisher distribution with central direction and apperture parameter kappa 
(Distribition size higher kappa = less variation)

-> Position

Type 0 - Uniform square with according to the size

Type 1 - Gaussian in all 3 directions with standard deviation specified by size

-> Momentum

Type 0 - Uniform square with according to the size

Type 1 - Gaussian in all 3 directions with standard deviation specified by size

Type 2 - Specific direction