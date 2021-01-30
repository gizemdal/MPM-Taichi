import taichi as ti
import numpy as np

ti.init(arch=ti.gpu) # Try to run on GPU

bunny_nodes = 0 # Number of nodes in the bunny OBJ mesh
n_particles, n_grid = 70000 * 2, 128
dx, inv_dx = 1 / n_grid, float(n_grid)
dt = 1e-4
p_vol, p_rho = (dx * 0.5)**2, 1
p_mass = p_vol * p_rho
E, nu = 1e3, 0.3 # Young's modulus and Poisson's ratio
mu_0, lambda_0 = E / (2 * (1 + nu)), E * nu / ((1+nu) * (1 - 2 * nu)) # Lame parameters - may change these later to model other materials

x = ti.Vector.field(3, dtype=float, shape=n_particles) # position
host_x = ti.Vector.field(3, dtype=float, shape=n_particles)
x_2d = ti.Vector.field(2, dtype=float, shape=n_particles) # 2d positions - this is necessary for circle visualization
v = ti.Vector.field(3, dtype=float, shape=n_particles) # velocity
C = ti.Matrix.field(3, 3, dtype=float, shape=n_particles) # affine velocity field
F = ti.Matrix.field(3, 3, dtype=float, shape=n_particles) # deformation gradient
material = ti.field(dtype=int, shape=n_particles) # material id
grid_v = ti.Vector.field(3, dtype=float, shape=(n_grid, n_grid, n_grid)) # grid node momentum/velocity
grid_m = ti.field(dtype=float, shape=(n_grid, n_grid, n_grid)) # grid node mass
gravity = ti.Vector.field(3, dtype=float, shape=()) # gravity

@ti.func 
def kirchoff_FCR(F, R, J, mu, la):
  return 2 * mu * (F - R) @ F.transpose() + ti.Matrix.identity(float, 3) * la * J * (J - 1) #compute kirchoff stress for FCR model (remember tau = P F^T)

@ti.kernel
def substep():
  # Step 1: clean grid data by zeroing out everything
  for i, j, k in grid_m:
    grid_v[i, j, k] = [0, 0, 0]
    grid_m[i, j, k] = 0

  # Particle state update and scatter to grid (P2G)
  for p in x: 
    
    # First for particle p, compute base index
    base = (x[p] * inv_dx - 0.5).cast(int)
    
    # Quadratic kernels  [http://mpm.graphics   Eqn. 123, with x=fx, fx-1,fx-2]
    fx = x[p] * inv_dx - base.cast(float)
    w = [0.5 * (1.5 - fx) ** 2, 0.75 - (fx - 1) ** 2, 0.5 * (fx - 0.5) ** 2]
    dw = [fx - 1.5, -2.0 * (fx - 1), fx - 0.5]

    mu, la = mu_0, lambda_0 #opportunity here to modify these to model other materials

    U, sig, V = ti.svd(F[p])
    J = 1.0

    for d in ti.static(range(3)):
      J *= sig[d, d]
    
    #Compute Kirchoff Stress
    kirchoff = kirchoff_FCR(F[p], U@V.transpose(), J, mu, la)

    #P2G for velocity and mass AND Force Update!
    for i, j, k in ti.static(ti.ndrange(3, 3, 3)): # Loop over 3x3x3 grid node neighborhood
      offset = ti.Vector([i, j, k])
      dpos = (offset.cast(float) - fx) * dx
      weight = w[i][0] * w[j][1] * w[k][2]
      
      # Compute for 3D
      dweight = ti.Vector.zero(float,3)
      dweight[0] = inv_dx * dw[i][0] * w[j][1] * w[k][2]
      dweight[1] = inv_dx * w[i][0] * dw[j][1] * w[k][2]
      dweight[2] = inv_dx * w[i][0] * w[j][1] * dw[k][2]
      
      force = -p_vol * kirchoff @ dweight # This is doing Step 6: Add elastic force

      # Step 2 & 3: Transfer mass and momentum from particles to grid
      grid_v[base + offset] += p_mass * weight * (v[p] + C[p] @ dpos) #momentum transfer
      grid_m[base + offset] += weight * p_mass #mass transfer

      grid_v[base + offset] += dt * force #add force to update velocity, don't divide by mass bc this is actually updating MOMENTUM
  
  # Gravity and Boundary Collision
  for i, j, k in grid_m:
    if grid_m[i, j, k] > 0: # No need for epsilon here
      # Step 4: Set velocity from momentum if mass != 0
      grid_v[i, j, k] = (1 / grid_m[i, j, k]) * grid_v[i, j, k] # Momentum to velocity
      # Step 5: Apply gravity on grid
      grid_v[i, j, k] += dt * gravity[None] * 9.8 # gravity
      
      #wall collisions - handle all 3 dimensions
      if i < 3 and grid_v[i, j, k][0] < 0:          grid_v[i, j, k][0] = 0 # Boundary conditions
      if i > n_grid - 3 and grid_v[i, j, k][0] > 0: grid_v[i, j, k][0] = 0
      if j < 3 and grid_v[i, j, k][1] < 0:          grid_v[i, j, k][1] = 0
      if j > n_grid - 3 and grid_v[i, j, k][1] > 0: grid_v[i, j, k][1] = 0
      if k < 3 and grid_v[i, j, k][2] < 0:          grid_v[i, j, k][2] = 0
      if k > n_grid - 3 and grid_v[i, j, k][2] > 0: grid_v[i, j, k][2] = 0
  
  # grid to particle (G2P)
  for p in x: 
    base = (x[p] * inv_dx - 0.5).cast(int)

    fx = x[p] * inv_dx - base.cast(float)
    w = [0.5 * (1.5 - fx) ** 2, 0.75 - (fx - 1.0) ** 2, 0.5 * (fx - 0.5) ** 2]
    dw = [fx - 1.5, -2.0 * (fx - 1), fx - 0.5]

    new_v = ti.Vector.zero(float, 3)
    new_C = ti.Matrix.zero(float, 3, 3)
    new_F = ti.Matrix.zero(float, 3, 3)

    for i, j, k in ti.static(ti.ndrange(3, 3, 3)): # loop over 3x3x3 grid node neighborhood
      dpos = ti.Vector([i, j, k]).cast(float) - fx
      g_v = grid_v[base + ti.Vector([i, j, k])]
      weight = w[i][0] * w[j][1] * w[k][2]

      # Compute for 3D
      dweight = ti.Vector.zero(float,3)
      dweight[0] = inv_dx * dw[i][0] * w[j][1] * w[k][2]
      dweight[1] = inv_dx * w[i][0] * dw[j][1] * w[k][2]
      dweight[2] = inv_dx * w[i][0] * w[j][1] * dw[k][2]

      new_v += weight * g_v
      new_C += 4 * inv_dx * weight * g_v.outer_product(dpos)
      new_F += g_v.outer_product(dweight)
    # Step 7: Interpolate new velocity back to particles
    v[p], C[p] = new_v, new_C
    # Step 8: Move the particles
    x[p] += dt * v[p] # advection
    x_2d[p] = [x[p][0], x[p][1]] # update 2d positions
    F[p] = (ti.Matrix.identity(float, 3) + (dt * new_F)) @ F[p] #updateF (explicitMPM way)

@ti.kernel
def reset():
  group_size = n_particles // 1
  for i in range(n_particles):
    # This is currently creating 2 bunnies
    x[i] = [host_x[i][0], host_x[i][1], host_x[i][2]]
    if i < n_particles // 2:
      material[i] = 0
      v[i] = [0.65, -0.9, 0]
    else:
      material[i] = 1
      v[i] = [-0.65, 0.9, 0]
    x_2d[i] = [x[i][0], x[i][1]]
    F[i] = ti.Matrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    C[i] = ti.Matrix.zero(float, 3, 3)
  
print("[Hint] Use WSAD/arrow keys to control gravity. Press R to reset.")
gui = ti.GUI("Explicit MPM", res=768, background_color=0x112F41)

# Load the bunny mesh

f = open("bunny_point.obj", "r") # This has 72027 vertices
for line in f:
  # Tokenize the line, we only care about lines of 5 tokens - "v (space) (x_pos) (y_pos) (z_pos)"
  tokens = line.split(' ')
  if len(tokens) != 4:
     continue
  if tokens[0] == 'v':
    # Write the position
    # We only allow positive coordinates for initial positions
    # We need to transform this mesh - scale it down and then translate each node
    host_x[bunny_nodes][0] = float(tokens[1]) * 0.2 + 0.2
    host_x[bunny_nodes][1] = float(tokens[2]) * 0.2 + 0.6
    host_x[bunny_nodes][2] = float(tokens[3][:-1]) * 0.2 + 0.2
    bunny_nodes += 1

# Copy the second bunny
for i in range(bunny_nodes):
  host_x[bunny_nodes + i][0] = host_x[i][0] + 0.5
  host_x[bunny_nodes + i][1] = host_x[i][1] - 0.4
  host_x[bunny_nodes + i][2] = host_x[i][2]

reset() # Call reset for the first time
gravity[None] = [0, -1, 0] # set initial gravity direction to -y

frame_idx = 0 # frame obj index

for frame in range(1000):
  if gui.get_event(ti.GUI.PRESS):
    if gui.event.key == 'r': 
      reset()
      frame_idx = 0 # reset the frame index as well
    elif gui.event.key in [ti.GUI.ESCAPE, ti.GUI.EXIT]: break
  if gui.event is not None: gravity[None] = [0, -1, 0] # if had any event
  if gui.is_pressed(ti.GUI.LEFT,  'a'): gravity[None][0] = -1
  if gui.is_pressed(ti.GUI.RIGHT, 'd'): gravity[None][0] = 1
  if gui.is_pressed(ti.GUI.UP,    'w'): gravity[None][1] = 1
  if gui.is_pressed(ti.GUI.DOWN,  's'): gravity[None][1] = -1
  for s in range(int(2e-3 // dt)):
    substep()
  colors = np.array([0xED553B,0x068587,0xEEEEF0], dtype=np.uint32)
  gui.circles(x_2d.to_numpy(), radius=1.5, color=colors[material.to_numpy()])
  # Write the positions into OBJ file
  #f = open('frames/frame' + str(frame_idx) + '.obj', 'w')
  #for p in x.to_numpy():
    #f.write("v %f %f %f\n" % (p[0], p[1], p[2]))
  #f.close()
  frame_idx += 1
  gui.show() # Change to gui.show(f'{frame:06d}.png') to write images to disk