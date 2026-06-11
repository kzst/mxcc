# mxcc

mxcc is an R package for Maxwell Control Charts in statistical process control. It provides functions to compute control limits and coefficients, calculate performance metrics (such as ARL, SDRL, MRL), and plot control charts specifically designed for monitoring Maxwell-distributed quality characteristics. The package is particularly useful in reliability engineering, physics, and other fields where data follows the Maxwell distribution.
---

## Key Features

* Computation of control limits for Maxwell Control Charts (V-chart, VSQ-chart, etc.).
* Calculation of coefficients and constants required for constructing Maxwell control charts.
* Performance evaluation using Average Run Length (ARL) and other metrics.
* Support for monitoring processes with Maxwell-distributed quality characteristics.
* Built-in plotting functions to visualize Maxwell control charts.
* Real Datasets are included in the package:
  

---

#### Author

* Zahid Khan
* Zsolt T. Kosztyan

#### Maintainer

* Zsolt T. Kosztyan

## Installation

Install the released version from CRAN :

install.packages("mxcc")


Or install the development version from GitHub:

```
library(devtools)
install_github("kzst/mxcc")
library(mxcc)

```


## Acknowledgement

This work has been implemented by the TKP2021-NVA-10 project with the support provided by the Ministry of Culture and Innovation of Hungary from the National Research, Development and Innovation Fund, financed under the 2021 Thematic Excellence Programme funding scheme.

