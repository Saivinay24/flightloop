# ALFA aircraft facts (approximate, with sources)

The CMU ALFA flights (Keipour, Mousaei, Scherer, IJRR 2021; doi 10.1184/R1/12707963) use a modified E-flite Carbon-Z T-28 with ArduPlane on a Pixhawk and an onboard computer. FlightLoop did not fly these aircraft; all numbers below are approximate.

| Quantity | Value used | Uncertainty | Source |
|---|---|---|---|
| Wingspan b | 1.98 m | +/-0.02 m | retailer spec tables, e.g. https://www.skyraccoon.com/aircraft/E-flite_Carbon-Z-T-28-Trojan_EFL013575 (77.95 in) |
| Wing area S | 0.692 m^2 | +/-5% | 1073 sq in listed by retailers (https://www.air-rc.com/aircraft/E-flite_Carbon-Z-T-28_EFL1350) |
| Mass m | 5.4 kg | +/-15% (ALFA adds an onboard computer and sensors; no mass is published) | Model Aviation review ready-to-fly weight 12 lb (https://www.modelaviation.com/node/1112) |
| Aspect ratio | 5.67 | from b and S | derived |
| Geometry for AVL | tapered wing c_root 0.437 m, c_tip 0.262 m, NACA 2412 stand-in, tail 0.104 m^2 at 0.85 m | approximate, not measured | `data/alfa/avl/t28_approx.avl` |

CD0 and k scale with m/S: a 15% mass error moves CL by 15% and the identified CD0 by up to 15%; this is propagated into the reported intervals (`results/sysid/alfa_summary.json`).
