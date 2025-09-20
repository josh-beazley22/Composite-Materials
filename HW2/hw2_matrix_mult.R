# --- Constants (ksi units) ---
E1 <- 18400    # ksi (18.4 Msi)
E2 <- 1600     # ksi (1.6 Msi)
G12 <- 950     # ksi (0.95 Msi)
v12 <- 0.28

sx <- 40       # ksi
sy <- -20      # ksi
sxy <- 20      # ksi

theta <- 30 * pi / 180  # convert degrees to radians

# --- Shorthand ---
c <- cos(theta)
s <- sin(theta)
Delta <- 1 - v12^2 * E2 / E1

# --- Matrices ---
R <- matrix(c(
  1, 0, 0,
  0, 1, 0,
  0, 0, 1/2
), nrow=3, byrow=TRUE)

Q <- matrix(c(
  E1/Delta, v12*E1/Delta, 0,
  v12*E1/Delta, E2/Delta, 0,
  0, 0, G12
), nrow=3, byrow=TRUE)

S <- matrix(c(
  1/E1, -v12/E1, 0,
  -v12/E1, 1/E2, 0,
  0, 0, 1/G12
), nrow=3, byrow=TRUE)

T <- matrix(c(
  c^2, s^2, 2*c*s,
  s^2, c^2, -2*c*s,
  -c*s, c*s, c^2 - s^2
), nrow=3, byrow=TRUE)

sigma <- matrix(c(sx, sy, sxy), nrow=3, byrow=TRUE)

# --- Compute strains ---
eps <- R %*% S %*% T %*% sigma

# --- Output ---
print("ε1, ε2, γ12 =")
print(eps)
