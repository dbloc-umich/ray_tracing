import numpy as np
import matplotlib.pyplot as plt

theta = np.linspace(0, 2*np.pi)
x = np.cos(theta)
y = np.sin(theta)
L = 0.5

# Case 1
plt.figure(figsize=(4.8, 2.4))
plt.plot(x,y)
plt.gca().annotate(r"$1$", xytext=(-1.5, 0), xy=(-1, 0), arrowprops=dict(arrowstyle="->"), va='center', fontsize=15)
plt.gca().annotate(r"$e^{-2\mu r}$", xytext=(1.5, 0), xy=(1, 0), arrowprops=dict(arrowstyle="<-"), va='center', fontsize=15)
plt.gca().set_aspect("equal")
plt.axis(False)
plt.xlim((-2, 2))
plt.ylim((-1, 1))
plt.savefig("../Figures/attenuation_case1.png")


# Case 2
nrows = 3
ncols = 2

fig, ax = plt.subplots(nrows, ncols)
fig.set_figheight(6.4*(nrows-0.2))
fig.set_figwidth(6.4*ncols)
fig.subplots_adjust(hspace=0.5, wspace=0.5)

for r in range(nrows):
    for c in range(ncols):
        ax[r][c].plot(x,y)
        ax[r][c].axis(False)
        ax[r][c].set_title(r*ncols+c+1)
        ax[r][c].set_aspect("equal")

# Initial ray
ax[0][0].annotate(r"$1$", xy=(-1, 0), xytext=(-1-L, 0), arrowprops=dict(arrowstyle="->"), va='center', fontsize=15)

# First step
ax[0][1].annotate(r"$T$", xytext=(-1+L, 0), xy=(-1, 0), arrowprops=dict(arrowstyle="<-"), va='center', fontsize=15)
ax[0][1].annotate(r"$R$", xytext=(-1-L, 0), xy=(-1, 0), arrowprops=dict(arrowstyle="<-", color="red"), va='center', color="red", fontsize=15)
ax[0][1].annotate(r"$TA$", xytext=(1-L, 0), xy=(1, 0), arrowprops=dict(arrowstyle="->"), va='center', fontsize=15)

ax[1][0].annotate(r"$T^2A$", xytext=(1+L/4, 0), xy=(1, 0), arrowprops=dict(arrowstyle="<-", color="red"), va='center', color="red", fontsize=15)
ax[1][0].annotate(r"$TRA$", xytext=(1-L, 0), xy=(1, 0), arrowprops=dict(arrowstyle="<-"), va='center', fontsize=15)
ax[1][0].annotate(r"$TRA^2$", xytext=(-1+L, 0), xy=(-1, 0), arrowprops=dict(arrowstyle="->"), va='center', fontsize=15)

ax[1][1].annotate(r"$TR^2A^2$", xytext=(-1+L, 0), xy=(-1, 0), arrowprops=dict(arrowstyle="<-"), va='center', fontsize=15)
ax[1][1].annotate(r"$T^2RA^2$", xytext=(-1-1.25*L, 0), xy=(-1, 0), arrowprops=dict(arrowstyle="<-", color="red"), va='center', color="red", fontsize=15)
ax[1][1].annotate(r"$TR^2A^3$", xytext=(1-1.25*L, 0), xy=(1, 0), arrowprops=dict(arrowstyle="->"), va='center', fontsize=15)

ax[2][0].annotate(r"$T^2R^2A^3$", xytext=(1+L/4, 0), xy=(1, 0), arrowprops=dict(arrowstyle="<-", color="red"), va='center', color="red", fontsize=15)
ax[2][0].annotate(r"$TR^3A^3$", xytext=(1-1.25*L, 0), xy=(1, 0), arrowprops=dict(arrowstyle="<-"), va='center', fontsize=15)
ax[2][0].annotate(r"$TR^3A^4$", xytext=(-1+L, 0), xy=(-1, 0), arrowprops=dict(arrowstyle="->"), va='center', fontsize=15)

ax[2][1].annotate(r"$TR^4A^4$", xytext=(-1+L, 0), xy=(-1, 0), arrowprops=dict(arrowstyle="<-"), va='center', fontsize=15)
ax[2][1].annotate(r"$T^2R^3A^4$", xytext=(-1-1.25*L, 0), xy=(-1, 0), arrowprops=dict(arrowstyle="<-", color="red"), va='center', color="red", fontsize=15)

plt.savefig("../Figures/attenuation_case2.png")

# Case 3
plt.figure(figsize=(4.8, 2.4))
plt.plot(x,y)

y0 = 0.5
x0 = np.sqrt(1-y0*y0)
plt.plot((-x0, x0), (y0, y0), color="black", linestyle="dashed")
plt.plot((0, 0), (0, y0), color="black", linestyle="dashed")
plt.text(0.1, y0/2, r"r/2")
plt.gca().annotate(r"$1$", xytext=(-x0-L, y0), xy=(-x0, y0), arrowprops=dict(arrowstyle="->"), va='center', fontsize=15)
plt.gca().annotate(r"$e^{-\sqrt{3}\mu r}$", xytext=(x0+L, y0), xy=(x0, y0), arrowprops=dict(arrowstyle="<-"), va='center', fontsize=15)
plt.gca().set_aspect("equal")
plt.axis(False)
plt.xlim((-2, 2))
plt.ylim((-1, 1))
plt.savefig("../Figures/attenuation_case3.png")