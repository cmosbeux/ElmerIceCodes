# Fluidity and Viscosity in Elmer

Let's define the viscosity as a function of the fluidity:

$$
\eta = \frac{1}{2} \left( EA \right)^{-1/3} 
$$

In Elmer, the fluidity uses the following unit:

$$ 
A_{\text{elmer}} = \left[ \text{MPa}^{-3} \cdot s^{-1} \right] 
$$

To convert $A_{\text{ref}}$ to this unit, we use:

$$
A_{\text{elmer}} = A_{\text{ref}} \left[ \text{Pa}^{-3} \cdot s^{-1} \right]  = A_{\text{ref}} \cdot 10^{18} \cdot \left( 31536 \cdot 10^{3} \right)  \left[ \text{MPa}^{-3} \cdot \text{a}^{-1} \right]
$$

If we consider $A_{\text{ref}} = 25 \cdot 10^{-25} \left[ \text{Pa}^{-3} \cdot s^{-1} \right]$, we get:

$$
A \approx 80 \left[ \text{MPa}^{-3} \cdot \text{s}^{-1} \right]
$$

We can then convert this to $\eta_{\text{elmer}}$:

$$
\eta_{\text{elmer}} = \frac{1}{2} \cdot A_{\text{elmer}}^{-1/3} = \frac{1}{2} \cdot 80^{-1/3} \approx 0.116 \left[ \text{MPa}^{3} \cdot \text{a}^{-1} \right]
$$

# Effect of viscosity and thickness change in the SSA framework

## 1D SSA Equation
In 1D, we can write the SSA equation as:
$$
\frac{\partial}{\partial x} \left( 4 H \nu \frac{\partial u}{\partial x} \right) = \rho g H \frac{\partial z_s}{\partial x}.
$$

## Viscosity Definition
The viscosity depends on the second invariant of the strain rate. In 1D, it only depends on the strain rate:
$$
\dot{\varepsilon} = \frac{\partial u}{\partial x}.
$$

The viscosity can be written as:
$$
\nu = \frac{1}{2} A^{-1/m} \left| \frac{\partial u}{\partial x} \right|^{\frac{1-n}{n}},
$$
where $m = \frac{1}{n} $.

## Viscosity Parameter
Let's define the viscosity parameter:
$$
B = A^{-1/n}.
$$

We can therefore write:
$$
\frac{\partial}{\partial x} \left( 4 H B \left| \frac{\partial u}{\partial x} \right|^{\frac{1-n}{n}} \frac{\partial u}{\partial x} \right) = \rho g H \frac{\partial z_s}{\partial x}.
$$

## Assumptions for Simplification
If we assume:
$$
\frac{\partial}{\partial x} \sim \frac{1}{L},
$$
we can rewrite the equation as:
$$
\frac{1}{L} \left( 4 H B \frac{u}{L} \left| \frac{u}{L} \right|^{\frac{1-n}{n}} \right) = \rho g H \frac{z_s}{L}.
$$

Simplifying further:
\[
\frac{H}{B} \left( \frac{1}{L} \right)^{1/m} \propto \rho g H \frac{H}{10}.
\]

## Scaling Relationship
From this, we infer:
$$
\mu \propto \left( \frac{H}{B} \right)^{1/n} L.
$$


