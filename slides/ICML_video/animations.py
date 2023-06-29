%matplotlib
import torch
import torch.nn as nn
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

class NN(nn.Module):
    def __init__(self):
        super().__init__()
        self.layers = nn.Sequential(
            nn.Linear(1, 256),
            nn.GELU(),
            *[ layer for _ in range(5) for layer in [nn.Linear(256, 256), nn.GELU() ]],
            nn.Linear(256, 1),
        )

    def forward(self, t):
        return self.layers(t.unsqueeze(1)).squeeze()


# True eqn:
# y'' + 2y' + 10y = 0
# y = e^(-t) * cos(3t)

y = lambda t: (-t).exp() * (3*t).cos()
full_t = torch.linspace(0, 5, 100)
full_y = y(full_t)
true_t = torch.linspace(0,5,5)
true_y = y(true_t) + torch.randn(5) * 0.05

fig, ax = plt.subplots(dpi = 300)
line, = ax.plot(full_t, full_y, linestyle = ":")
line.set_data([], [])
pts = ax.scatter([], [])

def update(frame):
    line.set_data(full_t[:frame], full_y[:frame])
    if frame < len(full_t):
        mask = true_t <= full_t[frame]
        cut_pts = list(zip(true_t[mask], true_y[mask]))
        pts.set_offsets(cut_pts)

ani = FuncAnimation(fig, update, frames = 150, interval = 20)
ani.save("assets/true_ani.gif")

## Naive nn ##
net = NN()
optim = torch.optim.Adam(net.parameters(), lr = 1e-3)

estim, = ax.plot(full_t, net(full_t).detach())
epoch = 0
ax.set_title(f"Epoch {epoch}")

def update(frame):
    global epoch
    for _ in range(1):
        out = net(true_t)
        loss = nn.functional.mse_loss(out, true_y)
        print(f"{loss=}")
        epoch = epoch + 1

        loss.backward()
        optim.step()
        optim.zero_grad()

    estim.set_data(full_t, net(full_t).detach())
    ax.set_title(f"Epoch {epoch}")

    return estim,

ani = FuncAnimation(fig, update, frames=100, interval = 100, repeat = False)
ani.save("assets/NN.gif")
plt.close(fig)







## PINN ##
fig, ax = plt.subplots(dpi=300)
ax.plot(full_t, full_y, linestyle = ":")
ax.scatter(true_t, true_y)

net = NN()
optim = torch.optim.Adam(net.parameters(), lr = 1e-3)

estim, = ax.plot(full_t, net(full_t).detach())
epoch = 0
ax.set_title(f"Epoch {epoch}")


def update(frame):
    global epoch
    for _ in range(5):
        out = net(true_t)
        loss_data = nn.functional.mse_loss(out, true_y)

        colloc = 5 * torch.rand(50, requires_grad = True)
        # colloc = true_t.clone().detach().requires_grad_(True)

        y = net(colloc)
        d, = torch.autograd.grad(y.sum(), colloc, create_graph = True)
        dd, = torch.autograd.grad(d.sum(), colloc, create_graph = True)

        loss_pde = (dd + 2*d + 10 * y).square().mean()

        loss = loss_data + loss_pde

        loss.backward()
        optim.step()
        optim.zero_grad()

        epoch += 1

    print(f'{loss_pde.detach()=} {loss_data.detach()=} {epoch=}')

    estim.set_data(full_t, net(full_t).detach())
    ax.set_title(f"Epoch {epoch}")

    return estim,

ani = FuncAnimation(fig, update, frames=100, interval = 100)
ani.save("assets/PINN.gif")
plt.close(fig)




## EPGP
fig, ax = plt.subplots(dpi=300)
ax.plot(full_t, full_y, linestyle = ":")
ax.scatter(true_t, true_y)

Y = true_y.unsqueeze(1)
Phi = torch.exp(-true_t) * torch.stack([ torch.cos(3 * true_t), torch.sin(3 * true_t) ])
Phi_star = torch.exp(-full_t) * torch.stack([ torch.cos(3 * full_t), torch.sin(3 * full_t) ])
log_s0 = torch.tensor(0., requires_grad = True)
log_s = torch.tensor([0., 0.], requires_grad = True)
n = len(true_t)
mr = 2

optim = torch.optim.Adam([log_s0, log_s], lr = 0.5)

Sigma = torch.diag_embed(log_s.exp())
Sigma_inv = torch.diag_embed((-log_s).exp())
A = mr * log_s0.exp() * Sigma_inv + Phi @ Phi.T
pred = (Phi_star.T @ A.inverse() @ Phi @ Y).squeeze().detach()

estim, = ax.plot(full_t, pred)
epoch = 0
ax.set_title(f"Epoch {epoch}")

def update(frame):
    global epoch
    for _ in range(1):
        print(f"{log_s0=} {log_s=}")
        Sigma = torch.diag_embed(log_s.exp())
        Sigma_inv = torch.diag_embed((-log_s).exp())
        A = mr * log_s0.exp() * Sigma_inv + Phi @ Phi.T
        print(f'{A=}')
        pred = (Phi_star.T @ A.inverse() @ Phi @ Y).squeeze().detach()

        loss = .5 * (-log_s0).exp() * (Y.T @ Y - Y.T @ Phi.T @ A.inverse() @ Phi @ Y) + (n - mr)/2 * log_s0 + .5 * log_s.sum() + .5 * A.det().log()
        print(f'{loss=}')

        loss.backward()
        optim.step()
        optim.zero_grad()

        epoch += 1

    estim.set_data(full_t, pred)
    ax.set_title(f"Epoch {epoch - 1}")

    return estim,

ani = FuncAnimation(fig, update, frames=20, interval = 100)
ani.save("assets/EPGP.gif")
plt.close(fig)








