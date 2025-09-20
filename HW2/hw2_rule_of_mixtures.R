
## Constants from Table 3.1 and 3.2
EF1 = 32.0 * 10^3
EF2 = 2.0 * 10^3
Em = 0.5 * 10^3
vF12 = 0.2
vm = 0.35
GF = 1.3 * 10^3
Gm = x

## Lookup values from Table 2.2
l.Vf = 0.65
l.E1 = 19.0 * 10^3
l.E2 = 1.5 * 10^3
l.G12 = 1.0 * 10^3
l.v12 = 0.22

ROM.E1 <- function(Vf) {
  ## parallel combination
  E1 = EF1 * Vf + Em * (1 - Vf)
  E1
}

ROM.v12 <- function(Vf) {
  ## parallel combination
  v12 = vF12 * Vf + vm * (1 - Vf)
  v12
}

ROM.E2 <- function(Vf) {
  ## series combination
  E2.inv = Vf / EF2 + (1-Vf) / Em
  1 / E2.inv
}

ROM.G12 <- function(Vf) {
  G12.inv = Vf / GF + (1-Vf) / Gm
  1 / G12.inv
}

RVE.E2 <- function(Vf) {
  ## 
  E2 = Em * (1 - sqrt(Vf) + sqrt(Vf) / (1 - sqrt(Vf) * (1 - Em/EF2)))
  E2
}

RVE.G12 <- function(Vf) {
  ## 
  G12 = Gm * (1 - sqrt(Vf) + sqrt(Vf) / (1 - sqrt(Vf) * (1 - Gm/GF)))
  G12
}

halpin.tsai.E2 <- function(Vf, xi=2) {
  eta = (EF2/Em - 1) / (EF2/Em + xi)
  E2 = Em * (1 + xi*eta*Vf) / (1 - eta*Vf)
  E2
}

halpin.tsai.G12 <- function(Vf, xi=1) {
  eta = (GF/Gm - 1) / (GF/Gm + xi)
  G12 = Gm * (1 + xi*eta*Vf) / (1 - eta*Vf)
  G12
}


par(mfrow = c(2,2))

# --- E1 ---
y.E1 <- c(ROM.E1(seq(0,1,len=100)), l.E1)
plot(NULL, xlim=c(0,1), ylim=range(y.E1, na.rm=TRUE),
     xlab="Fiber volume fraction", ylab="ksi", main="E1")
curve(ROM.E1(x), 0, 1, lty=1, add=TRUE)
points(l.Vf, l.E1, pch=16)
legend("bottomright", legend="Rule of Mixture", lty=1, bty="o", cex=0.8)

# --- v12 ---
y.v12 <- c(ROM.v12(seq(0,1,len=100)), l.v12)
plot(NULL, xlim=c(0,1), ylim=range(y.v12, na.rm=TRUE),
     xlab="Fiber volume fraction", ylab="unitless", main="v12")
curve(ROM.v12(x), 0, 1, lty=1, add=TRUE)
points(l.Vf, l.v12, pch=16)
legend("topright", legend="Rule of Mixture", lty=1, bty="o", cex=0.8)

# --- E2 ---
y.E2 <- c(ROM.E2(seq(0,1,len=100)),
          RVE.E2(seq(0,1,len=100)),
          halpin.tsai.E2(seq(0,1,len=100)),
          l.E2)
plot(NULL, xlim=c(0,1), ylim=range(y.E2, na.rm=TRUE),
     xlab="Fiber volume fraction", ylab="ksi", main="E2")
curve(ROM.E2(x), 0, 1, lty=1, add=TRUE)
curve(RVE.E2(x), 0, 1, lty=2, add=TRUE)
curve(halpin.tsai.E2(x), 0, 1, lty=3, add=TRUE)
points(l.Vf, l.E2, pch=16)
legend("topleft", legend=c("Rule of Mixture","RVE","Halpin-Tsai"),
       lty=1:3, bty="o", cex=0.8)

# --- G12 ---
y.G12 <- c(ROM.G12(seq(0,1,len=100)),
           RVE.G12(seq(0,1,len=100)),
           halpin.tsai.G12(seq(0,1,len=100)),
           l.G12)
plot(NULL, xlim=c(0,1), ylim=range(y.G12, na.rm=TRUE),
     xlab="Fiber volume fraction", ylab="ksi", main="G12")
curve(ROM.G12(x), 0, 1, lty=1, add=TRUE)
curve(RVE.G12(x), 0, 1, lty=2, add=TRUE)
curve(halpin.tsai.G12(x), 0, 1, lty=3, add=TRUE)
points(l.Vf, l.G12, pch=16)
legend("topleft", legend=c("Rule of Mixture","RVE","Halpin-Tsai"),
       lty=1:3, bty="o", cex=0.8)
