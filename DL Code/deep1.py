import torch
import torch.nn as nn
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import time

# Start the timer
start_time = time.time()

# SETUP
def load_csv_matrix(file_path):
    df = pd.read_csv(file_path, header=None, skiprows=1)
    matrix = df.values.astype(np.float32)
    return torch.tensor(matrix, dtype=torch.float32)

S_true = load_csv_matrix("Sm.csv")
V_true = load_csv_matrix("Vm.csv")
E_true = load_csv_matrix("Em.csv")
I_true = load_csv_matrix("Im.csv")
R_true = load_csv_matrix("Rm.csv")


x_csv = torch.arange(-3.0, 3.01, 0.1, dtype=torch.float32)  # 61
t_csv = torch.arange(0, 61, 6,dtype=torch.float32)        # 11
X_ref, T_ref = torch.meshgrid(x_csv, t_csv, indexing="ij")

def prepare_supervised_data(ref_tensor,name, num_samples=20):
    X_flat = X_ref.reshape(-1, 1)
    T_flat = T_ref.reshape(-1, 1)
    Val_flat = ref_tensor.reshape(-1, 1)

    XT_all = torch.cat([X_flat, T_flat], dim=1)
    indices = torch.randperm(XT_all.shape[0])[:num_samples]
    XT_sampled = XT_all[indices]
    Val_sampled = Val_flat[indices]

    return XT_sampled, Val_sampled

data_S, target_S = prepare_supervised_data(S_true, "S")
data_V, target_V = prepare_supervised_data(V_true, "V")
data_E, target_E = prepare_supervised_data(E_true, "E")
data_I, target_I = prepare_supervised_data(I_true, "I")
data_R, target_R = prepare_supervised_data(R_true, "R")

#____________________________________________________________________________________

# XT is a tensor of shape [N, 2] containing (x, t) pairs

x = torch.linspace(-3, 3 , 70 ,dtype=torch.float32)
t = torch.linspace(0, 60 , 70,dtype=torch.float32)
X,T = torch.meshgrid(x, t, indexing="ij")
XT = torch.stack([X.reshape(-1), T.reshape(-1)], dim=1).requires_grad_(True)


#Constants
b = torch.tensor(0.5140)
be=torch.tensor(0.25)
bi=torch.tensor(1)
bv=torch.tensor(0.9)
# rho = torch.tensor(0.1)
phi = torch.tensor(1/20)
sigma = torch.tensor(1/2)
gamma = torch.tensor(1/5)
delta=torch.tensor(1/365)
theta = torch.tensor(1/365)
# mu=torch.tensor(5.50*(10)**-8)
r=torch.tensor(7.140*(10)**-5)
k=torch.tensor(1.857*(10)**-4)
alpha=torch.tensor(9.300*(10)**-6)

d1 = torch.tensor(0.05)
d2 = torch.tensor(0.05)
d3 = torch.tensor(0.025)
d4 = torch.tensor(0.001)
d5 = torch.tensor(0.0)
#___________________________________________________________________________________


# PINN model
class PINN(nn.Module):
    def __init__(self):
        super().__init__()
        self.layer = nn.Sequential(
            nn.Linear(2, 100),
            nn.SiLU(),
            nn.Linear(100, 100),
            nn.SiLU(),
            nn.Linear(100, 100),
            nn.SiLU(),
            nn.Linear(100, 100),
            nn.SiLU(),
            nn.Linear(100, 5),
            nn.Softplus()
                            )


    def forward(self, x):
        out= self.layer(x)
        return out

class Net:
    def __init__(self):
        self.model = PINN()
        

        self.XT = XT.clone().requires_grad_(True)
        self.criterion = nn.MSELoss()
        self.ic_mask = (self.XT[:, 1] == 0.0)
        self.bc_mask = (self.XT[:, 0] == -3.0) | (self.XT[:, 0] == 3.0)

        self.data_S, self.target_S = data_S, target_S
        self.data_V, self.target_V = data_V, target_V
        self.data_E, self.target_E = data_E, target_E
        self.data_I, self.target_I = data_I, target_I
        self.data_R, self.target_R = data_R, target_R

        self.criterion = torch.nn.MSELoss()
        self.epoch = 1


        self.optimizer = torch.optim.LBFGS(
            self.model.parameters(),
            lr=1.0,
            max_iter=50000,
            max_eval=50000,
            history_size=50,
            tolerance_grad=1e-3,
            tolerance_change=1.0 * np.finfo(float).eps,
            line_search_fn="strong_wolfe",   # better numerical stability
        )

        self.adam = torch.optim.Adam(self.model.parameters())

    def compute_residuals(self):
        y = self.model(self.XT)


        if self.epoch == 1000:
            print("output row of model:\n", y[0])

        S, V, E, I, R = y[:, 0:1], y[:, 1:2], y[:, 2:3], y[:, 3:4], y[:, 4:5]
       
        grad_S = torch.autograd.grad(S, self.XT, grad_outputs=torch.ones_like(S), create_graph=True)[0]
        dS_dx, dS_dt = grad_S[:, 0:1], grad_S[:, 1:2]
        d2S_dx2 = torch.autograd.grad(dS_dx, self.XT, grad_outputs=torch.ones_like(dS_dx), create_graph=True)[0][:, 0:1]

        grad_V = torch.autograd.grad(V, self.XT, grad_outputs=torch.ones_like(V), create_graph=True)[0]
        dV_dx, dV_dt = grad_V[:, 0:1], grad_V[:, 1:2]
        d2V_dx2 = torch.autograd.grad(dV_dx, self.XT, grad_outputs=torch.ones_like(dV_dx), create_graph=True)[0][:, 0:1]

        grad_E = torch.autograd.grad(E, self.XT, grad_outputs=torch.ones_like(E), create_graph=True)[0]
        dE_dx, dE_dt = grad_E[:, 0:1], grad_E[:, 1:2]
        d2E_dx2 = torch.autograd.grad(dE_dx, self.XT, grad_outputs=torch.ones_like(dE_dx), create_graph=True)[0][:, 0:1]

        grad_I = torch.autograd.grad(I, self.XT, grad_outputs=torch.ones_like(I), create_graph=True)[0]
        dI_dx, dI_dt = grad_I[:, 0:1], grad_I[:, 1:2]
        d2I_dx2 = torch.autograd.grad(dI_dx, self.XT, grad_outputs=torch.ones_like(dI_dx), create_graph=True)[0][:, 0:1]

        grad_R = torch.autograd.grad(R, self.XT, grad_outputs=torch.ones_like(R), create_graph=True)[0]
        dR_dx, dR_dt = grad_R[:, 0:1], grad_R[:, 1:2]
        d2R_dx2 = torch.autograd.grad(dR_dx, self.XT, grad_outputs=torch.ones_like(dR_dx), create_graph=True)[0][:, 0:1]

        res_S = dS_dt + b * bi * S * I + b * be * S * E - alpha * I * S + phi * S + r * S - delta * R - d1 * d2S_dx2 - theta * V - r
        res_V = dV_dt + b * be * bv * E * V + b * bi * bv * I * V - alpha * I * V - r * V + theta * V - phi * S - d2 * d2V_dx2
        res_E = dE_dt - b * be * E * S - b * bi * I * S - b * be * bv * E * V - b * bi * bv * I * V - alpha * I * E + (r + k + sigma) * E - d3 * d2E_dx2
        res_I = dI_dt - sigma * E + (r + alpha + gamma) * I - alpha * I*2 - d4* d2I_dx2
        res_R = dR_dt - k * E - gamma * I + r * R + delta * R - alpha * I * R - d5 * d2R_dx2


        return res_S, res_V, res_E, res_I, res_R



    def loss_func(self):

        self.adam.zero_grad()
        self.optimizer.zero_grad()
        res = self.compute_residuals()
        loss_pde = sum(self.criterion(r, torch.zeros_like(r)) for r in res)

        x_ic = self.XT[self.ic_mask, 0:1]
        t_ic = self.XT[self.ic_mask, 1:2]
        XT_ic = torch.cat([x_ic, t_ic], dim=1)

        S0 = 0.86 * torch.exp(-((x_ic / 1.4) ** 2))
        V0 = 0.10 * torch.exp(-((x_ic / 1.4) ** 2))
        E0 = 0.00 * torch.exp(-((x_ic / 1.0) ** 2))
        I0 = 0.04 * torch.exp(-(x_ic ** 2))
        R0 = 0.00 * torch.exp(-(x_ic ** 2))
        Y0 = torch.cat([S0, V0, E0, I0, R0], dim=1)

        pred_ic = self.model(XT_ic)
        loss_ic = self.criterion(pred_ic, Y0)

        XT_bc = self.XT[self.bc_mask].detach().clone().requires_grad_(True)
        output_bc = self.model(XT_bc)
        d_bc = [torch.autograd.grad(output_bc[:, i:i+1], XT_bc, grad_outputs=torch.ones_like(output_bc[:, i:i+1]), create_graph=True)[0][:, 0:1] for i in range(5)]
        loss_bc = sum(self.criterion(deriv, torch.zeros_like(deriv)) for deriv in d_bc)


        pred_S = self.model(self.data_S)[:, 0:1]
        pred_V = self.model(self.data_V)[:, 1:2]
        pred_E = self.model(self.data_E)[:, 2:3]
        pred_I = self.model(self.data_I)[:, 3:4]
        pred_R = self.model(self.data_R)[:, 4:5]

        loss_data = (
            self.criterion(pred_S, self.target_S) +
            self.criterion(pred_V, self.target_V)*1.50 +
            self.criterion(pred_E, self.target_E) +
            self.criterion(pred_I, self.target_I) +
            self.criterion(pred_R, self.target_R)
        )


        loss =  loss_pde +  loss_ic +  loss_bc + 0.1*loss_data


        loss.backward()
        if self.epoch % 100 == 0:
            print(f"Epoch {self.epoch}, Loss: {loss.item():.6f}, data: {loss_data.item():.6f},PDE: {loss_pde.item():.6f}, IC: {loss_ic.item():.6f}, BC: {loss_bc.item():.6f}")
        self.epoch = self.epoch + 1
        return loss

    def train(self):
        self.model.train()
        for epoch in range(1000):
            self.adam.step(self.loss_func)
            self.adam.zero_grad()
        self.optimizer.step(self.loss_func)
        self.optimizer.zero_grad()
    def eval_(self):
        self.model.eval()


net = Net()
net.train()


# End the timer
end_time = time.time()
# Calculate the duration
elapsed_time = end_time - start_time

# === PLOTTING CODE ===
net.eval_()
with torch.no_grad():
    prediction = net.model(net.XT)
    S_pred = prediction[:, 0].reshape(len(x), len(t))
    V_pred = prediction[:, 1].reshape(len(x), len(t))
    E_pred = prediction[:, 2].reshape(len(x), len(t))
    I_pred = prediction[:, 3].reshape(len(x), len(t))
    R_pred = prediction[:, 4].reshape(len(x), len(t))


desired_times = np.arange(0, 61, 6)  # [0, 6, 12, ..., 60]
plot_time_indices = [int(np.argmin(np.abs(t.numpy() - dt))) for dt in desired_times]

colors = plt.cm.viridis(np.linspace(0, 1, len(plot_time_indices)))

def plot_quantity(quantity, label):
    plt.figure(figsize=(10, 6))
    for i_line, i_full_t in enumerate(plot_time_indices):
        plt.plot(x.numpy(), quantity[:, i_full_t].numpy(), color=colors[i_line], label=f"t={t[i_full_t].item():.0f}")
    plt.title(f"{label}(x, t) over X ")
    plt.xlabel("x")
    plt.ylabel(f"{label}(x, t)")
    plt.grid(True)
    plt.legend(loc="upper right", fontsize="small")
    plt.tight_layout()
    plt.show()

# Call plots
plot_quantity(S_pred, "S")
plot_quantity(V_pred, "V")
plot_quantity(E_pred, "E")
plot_quantity(I_pred, "I")
plot_quantity(R_pred, "R")


# === PLOTTING AT x = 0 OVER TIME — SEPARATE GRAPHS ===

# Find the index where x == 0
x0_index = torch.argmin(torch.abs(x)).item()

# Get values at x = 0 for each quantity
S_x0 = S_pred[x0_index, :].numpy()
V_x0 = V_pred[x0_index, :].numpy()
E_x0 = E_pred[x0_index, :].numpy()
I_x0 = I_pred[x0_index, :].numpy()
R_x0 = R_pred[x0_index, :].numpy()

# Plot each variable at x = 0 in a separate graph
variables = {
    "S(x=0,t)": S_x0,
    "V(x=0,t)": V_x0,
    "E(x=0,t)": E_x0,
    "I(x=0,t)": I_x0,
    "R(x=0,t)": R_x0
}

for label, values in variables.items():
    plt.figure(figsize=(8, 5))
    plt.plot(t.numpy(), values, label=label)
    plt.xlabel("Time (t)")
    plt.ylabel(label)
    plt.title(f"{label} over time at x = 0")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()


#error part
# from model
# Define x and t values
t_values = torch.arange(0, 61, 6)       # t = 0, 6, ..., 60 → shape = [11]
x_values = torch.arange(-3.0, 3.01, 0.1)  # x = -3.0, -2.9, ..., 3.0 → shape = [61]

# Find indices directly
t_indices = torch.tensor([torch.argmin(torch.abs(t - val)).item() for val in t_values])
x_indices = torch.tensor([torch.argmin(torch.abs(x - val)).item() for val in x_values])

S_pred_tensor = S_pred.index_select(0, x_indices).index_select(1, t_indices)
V_pred_tensor = V_pred.index_select(0, x_indices).index_select(1, t_indices)
E_pred_tensor = E_pred.index_select(0, x_indices).index_select(1, t_indices)
I_pred_tensor = I_pred.index_select(0, x_indices).index_select(1, t_indices)
R_pred_tensor = R_pred.index_select(0, x_indices).index_select(1, t_indices)


#time
print(f"Time taken: {elapsed_time:.4f} seconds")

# true data
# function
def load_csv_matrix(file_path):
    df = pd.read_csv(file_path, header=None, skiprows=1)
    matrix = df.values.astype(np.float32)
    return torch.tensor(matrix, dtype=torch.float32)

# data load
S_true = load_csv_matrix("Sm.csv")
V_true = load_csv_matrix("Vm.csv")
E_true = load_csv_matrix("Em.csv")
I_true = load_csv_matrix("Im.csv")
R_true = load_csv_matrix("Rm.csv")


#error calculation
# error matrix calculation
S_diff = S_pred_tensor - S_true
V_diff = V_pred_tensor - V_true
E_diff = E_pred_tensor - E_true
I_diff = I_pred_tensor - I_true
R_diff = R_pred_tensor - R_true

#error matrix print
print("S Error Matrix:\n", S_diff)
print("V Error Matrix:\n", V_diff)
print("E Error Matrix:\n", E_diff)
print("I Error Matrix:\n", I_diff)
print("R Error Matrix:\n", R_diff)

# Mean
mean_error_S = torch.mean(torch.abs(S_diff)).item()
mean_error_V = torch.mean(torch.abs(V_diff)).item()
mean_error_E = torch.mean(torch.abs(E_diff)).item()
mean_error_I = torch.mean(torch.abs(I_diff)).item()
mean_error_R = torch.mean(torch.abs(R_diff)).item()


#print mean
print(f"\nMean Error for S: {mean_error_S:.6f}")
print(f"Mean Error for V: {mean_error_V:.6f}")
print(f"Mean Error for E: {mean_error_E:.6f}")
print(f"Mean Error for I: {mean_error_I:.6f}")
print(f"Mean Error for R: {mean_error_R:.6f}")
