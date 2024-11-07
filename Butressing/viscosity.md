## Fluidity and Viscosity in Elmer

Let's define the viscosity as a function of the fluidity:

$$
\eta = \frac{1}{2} \left( EA \right)^{1/3} 
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
