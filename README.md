# zjuchenGeochem

This package contains several functions to facilitate your geochemical research.
The functions will be irregularly implemented.
Currently the function list includes:

### `seawater_density`:

**Function**: Calculate the seawater density at given P-T-X condition.

**Usage**: `seawater_density (Temperature = 100, P = 20, S = 35)`

**Source**: (P,ρ, T) properties of standard seawater.
by Safarov et al., (2012, 2013)

**Ranges**: absolute salinity SA=0-55.529 g kg^-1^, *T* = 273.15 - 468.15 K, *P* up to 140 MPa (1400 bar).

**Equations**:

$$P=A(\rho)^2+B(\rho)^8+C(\rho)^{12}$$ *ρ* in g cm^-3^, *P* in MPa, *T* in K

$$A=\Sigma^5_{i=1}*T^i*\Sigma^2_{j=0}*a_{ij}*S_A^j$$ $$B=\Sigma^4_{i=0}*T^i*\Sigma^2_{j=0}*b_{ij}*S_A^j$$ $$C=\Sigma^4_{i=0}*T^i*\Sigma^2_{j=0}*c_{ij}*S_A^j$$

The coefficients for a~ij~, b~ij~, and c~ij~ are shown in the Table 2 of Safarov et al., 2012.

### `critical_pressure`:

**Function**: Calculate the critical pressure of seawater at given temperature.

**Usage**: `critical_pressure(Temperature = 25, type = "saturated")`

`type = "saturated"`:

`type = "full"`:

**Source**: Critical curve of seawater proposed by Driesner and Heinrich, 2007.

**Equation 5a** for temperatures below critical temperature of pure water (373.976 ^o^C).

This equation is used for **pure water**, however, at temperature of below critical point of pure water, the critical curves of pure water and seawater are almost identical (Fig. 1)

$$P_{crit}=P^{H2O}_{crit}+\Sigma^7_{n=1}c_n(T^{H2O}_{crit}-T)^{c_nA}$$ *P~crit~* = 220.54915 bar.
The values for the constants are shown in Table 4.

### `quartz_sol`:

**Function**: Calculate the SiO~2~ concentration of the fluid equilibrates with quartz at given temperature.
The corresponding pressure at each temperature is calculated using the `critical_pressure` function.

**Usage**: `quartz_sol(Temperature = 100, S = 35, ref = "vonDamm1991")`

`ref = "vonDamm1991"`.
Equation derived from Von Damm et al., (1991).

$$lnM=a+bln\rho+(C+dT^2)/T+eP/T$$where *M* is the molarity of dissolved SiO~2~, *ρ* is the density of the solution (g cm^-3^), *P* is the pressure (bar), *T* is temperature in *K*.

a=-2.32888, b = 1.79547, c=-2263.62, d=0.00407350, e=0.0398808.

This equation is proposed to be valid at temperature range from 45 ^o^C to 900 ^o^C and pressure up to 9860 bars.

`ref = "Gunnarsson2000"`.
Equation derived from Gunnarsson and Arnorsson (2000).

$$logK_{am.silica}=-8.476-485.24 \times T^{-1} - 2.268 \times 10^{-6} \times T^2+3.068\times logT$$ $$logK_{quartz}=-34.188+197.47 \times T^{-1} - 5.851 \times 10^{-6} \times T^2+12.245\times logT$$ *T* in K and valid at 0 - 350 ^o^C at 1 bar below 100 ^o^C and P~sat~ at higher temperatures.
The calculated SiO~2~ concentration is in an unit of ppm, and late transformed into mM.

**The Von Damm's equation is preferably used**

### `dissolution_constant`:

**Function**: Calculate the dissolution constant of a specified mineral at given temperature

**Usage**: `dissolution_constant (Temperature = 100, mineral = "anhydrite")`

`mineral = "anhydrite"`: Anhydrite dissolution by Arnorsson et al., (1982).
Table 5

$$CaSO_4(s)=Ca^{2+}+SO_4^{2-}$$ $$log10(K)=3.94-677/T-12.39*T/1000$$ This equation is valid at 0-370 ^o^C.

`mineral = "calcite"`: Table 5 in Arnorsson et al., (1982).
The temperature dependence of **apparent** equilibrium constants of calcite is:

$$logK_{calcite}=log(m_{Ca^{2+}} \times
m_{CO_3^{2-}})=-1.46-41/T-17.41\times 10^{-6}T^2$$

where *T* is absolute temperature.

`mineral = "barite"`

`mineral = "celestite"`

The solubilities of barite and celestite are given as the following equation [\@monnin1999]:

$$MSO_4(s) \rightleftharpoons M^{2+}+SO_4^{2-}$$

$$K_{sp}=m_{M^{2+}}\times m_{SO_4^{2-}}\times\gamma_{M^{2+}}\times\gamma_{SO_4^{2-}}$$ where *m* is the molality and $\gamma$ the activity coefficient of the designated aqueous species.

Activity coefficients are generated from classic Debye-Huckel theory using the *R* package *CHNOSZ* developed by [@dick2019].

The temperature dependence of the dissolution constants for barite and celestite at saturation pressure is [@monnin1999]:

$$lnK_{sp,barite}=275.053-43.014lnT-15806.30/T$$

$$lnK_{sp, celestite}=224.069-35.9422lnT-10302.32/T$$

The uncertainty is 1% at 298.15K and of 10% at other temperatures for the barite dissolution.

## `Henry_constant`

**Function**: Calculate the Henry's law constant of a specified gas at given temperature.
There are several values to express Henry's law constant.

The Henry volatility defined via aqueous-phase mixing ratio: $K^{px}_H=\frac {p}{x}$ with a unit of atm.
The higher of this value, the lower proportion of gas dissolution.

The dimensionless Henry solubility: $H^{cc}=\frac{c_{aq}}{c_{gas}}$ which is dimensionless.
The higher of this value, the higher gas dissolution degree.

This function will return the $K^{px}_H$ value of a gas at given temperature.

**Usage**: `Henry_constant (Temperature = 25, gas = "Ar")`

`gas = "Ar"`: According to Krause and Benson (1989).
$$T^{*2}ln(K_H)=T^{*2}A_0+A_1(1-T^*)^{1/3}+A_2(1-T^*)^{2/3}$$ $T^*\equiv T/T_{cl}$, *T~cl~* is the critical point of water (647 K).
*A~0~*=4.289125, *A~1~*=19.225988, *A~2~*=-21.603721.

**Warning: This equation is derived from experiments at 0 - 60 ^o^C.** But t fits the high temperature data by other researchers.

`gas = "N2"`:

$$K^{px}_H=\frac{1}{10^{-pB}/55.55}$$ $$pB=9.755-\frac{1121}{T}-0.009169*T$$ where T is temperature in K.
**This equation is valid at 0 -350 ^o^C**.

`gas = "N2/Ar"`: return the molar ratio of N~2~/Ar dissolved in the fluid that equilibrate with air.

$$Ar_{aq}\sim Ar_{air}/K_{Ar}$$ $$N2_{aq}\sim N2_{air}/K_{N2}$$ $$(\frac{N_2}{Ar})_{aq}=(\frac{N_2}{Ar})_{air}*(\frac{K_{Ar}}{K_{N2}})$$

### `Horita, 2001`:

**Function**: Calculate the equilibration temperature of CO~2~ and CH~4~ using the equation proposed by Horita (2001).

$$
10^3ln\alpha(CO_2-CH_4)=0.16+11.754(10^6/T^2)-2.3655(10^9/T^3)+0.2054(10^{12}/T^4)
$$

*T* = 273 - 1573 K, 1$\sigma$ = ± 0.21‰, n = 44

**Usage**: `Horita_T (Temperature = 25)`, input temperature, return $\Delta_{CO_2-CH_4}$

`Horita (dC=25)`: input $\Delta_{CO_2-CH_4}$, return equilibration temperature

## `d18O_temperature`:

## `REE_normalization`:

### `CO2_dissolution`:

**Function**: Calculate the molar ratio between dissolved CO~2~ and that remained in the gas phase.

According to Henry's law (Sander, 2015), the dissolution of gas in an ideal state can be expressed as: $$H^{cc} = C_{aq}/C_{gas}$$

where *C~aq~* and *C~gas~* are the gas dissolved in the aqueous phase and remained in the gas phase, respectively.
*H^cc^* is the dimensionless Henry's law constants, which changes with temperature: $$H^{cc}=H_{\theta}^{cc}\times exp[\frac{-\Delta_{sol}H}{R}(\frac{1}{T}-\frac{1}{T_\theta})]$$ where $H^{cc}_{\theta}$ is the Henry's law constant at *T*$_\theta$ (298.15 K), *T* is the temperature of the dissolution system.
$-\Delta_{sol}H/R$ is a constant refers to the enthalpy of dissolution which does not change much with temperature.
The $H^{cc}_{\theta}$ and $-\Delta_{sol}H/R$ of CO~2~ are 0.83 and 2400 K, respectively.

Henry's law can be used for modeling the solubility of CO~2~ in aqueous solution for pressures up to 100 MPa (Carrol and Mather, 1992).
Therefore, the pressure posed little influence on the calculation.

On the other hand,the fraction of the total dissolved as CO~2~(aq) changes with pH, while the *C~aq~* is controlled by the Henry's law.

CO~2~(aq) + H~2~O ⇋ H~2~CO~3~

H~2~CO~3~ + H~2~O ⇋ H~3~O^+^ + HCO~3~^−^

*K~a1~ = 4.45 10^-7^* at 25 ^o^C

HCO~3~^−^ + H~2~O ⇋ H~3~O^+^ + CO~2~^−^

*K~a2~ = 4.67 10^-11^* at 25 ^o^C

Therefore, the total CO~2~(aq)/CO~2~(aq) = 1 + *K~a1~*/[H^+^] + *K~a1~*\*K~a2~/[H^+^]^2^

The *K~a~* values are temperature dependent.
According to Van't Hoff equation:

$$\frac{dln(K_a)}{dT} = \frac{\Delta H^{\theta}}{RT^2}$$ $$ln\frac{K_2}{K_1} = \frac{-\Delta H^{\theta}}{R} (\frac{1}{T_2} - \frac{1}{T_1})$$

Based on these theories, we could calculate the ratios between CO~2~ in aqueous solution and in gas phase.

**Usage**: `CO2_dissolution(Temperature = 100, pH = 8)`

The output is: 1.61824, indicating that the molar ratio between dissolved CO~2~ and gas phase CO~2~ is 1.62

## My favorite ggplot themes

```{r}
theme.chen2 <- theme(axis.title = element_text(size=18),
                     axis.text = element_text(size=14),
                     legend.text = element_text(size=12),
                     legend.title = element_text(size=14))

theme.chen1 <-  theme(
  axis.title = element_text(size = 18),
  axis.text = element_text(size = 14),
  legend.text = element_text(size = 12),
  legend.title = element_text(size = 14),
  axis.ticks.length = unit(0.15,"cm"),
  # axis.text.x.bottom = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")),
  # axis.text.y.left = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")),
  # axis.text.x.top = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")),
  # axis.text.y.right = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")),
  axis.line = element_line(size = 1),
  # axis.title.x.bottom = element_text(margin=unit(c(-.4,-.1,-.1,-.1), "cm")),
  # axis.title.y.left = element_text(margin=unit(c(-.1,-.4,-.1,-.1), "cm")),
  panel.background = element_blank()
)

insert_minor <- 
  function(major_labs, n_minor,n_breaks) {
    labs <- c( sapply( major_labs, 
                       function(x) c(x, rep("", n_breaks) ) ) )
                              labs[1:(length(labs)-n_minor)]
                              }
```
