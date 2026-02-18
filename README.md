# TWL

An R package for calculating the **Thermal Work Limit (TWL)** — the maximum sustainable metabolic rate (W/m²) for workers in thermally stressful environments. Implements the Brake & Bates (2002) methodology.

## Installation

Install from GitHub using the `remotes` or `pak` package:

```r
# Using remotes
install.packages("remotes")
remotes::install_github("vorpalvorpal/TWL")

# Using pak
install.packages("pak")
pak::pkg_install("vorpalvorpal/TWL")
```

Or install from a local clone:

```r
# Clone the repository, then:
install.packages(".", repos = NULL, type = "source")
```

### Dependencies

The package requires R (>= 4.1) and the following packages, which are installed automatically:

- `dplyr`
- `httr2`
- `cli`
- `purrr`

## Usage

### Basic calculation (all parameters supplied)

```r
library(TWL)

twl <- generate_twl(
  datetime     = as.POSIXct("2024-01-15 10:00:00", tz = "UTC"),
  latitude     = -31.95,
  longitude    = 115.86,
  temp         = 36,    # °C dry bulb
  wind_speed   = 0.5,   # m/s at body level
  RH           = 40,    # % relative humidity
  direct_solar = 600,   # W/m²
  diffuse_solar = 100,  # W/m²
  pressure     = 1013   # hPa
)
```

### Auto-fetch missing weather from Open-Meteo

Any weather parameter can be omitted and will be fetched automatically from the [Open-Meteo API](https://open-meteo.com/) (free, no API key required). An internet connection is required.

```r
twl <- generate_twl(
  datetime  = as.POSIXct("2024-01-15 10:00:00", tz = "UTC"),
  latitude  = -31.95,
  longitude = 115.86,
  temp      = 36
  # wind_speed, RH, direct_solar, diffuse_solar, pressure fetched from API
)
```

### Categorise results

```r
categorise_twl(twl)
# "Unrestricted", "Acclimatisation", "Buffer", or "Withdrawal"

twl_colour(twl)
# Hex colour codes for visualisation
```

## TWL Categories

| TWL (W/m²) | Category | Colour |
|---|---|---|
| ≥ 220 | Unrestricted | Green |
| 140 – 220 | Acclimatisation | Amber |
| 115 – 140 | Buffer | Orange |
| < 115 | Withdrawal | Red |

## Reference

Brake, D.J. and Bates, G.P. (2002) Limiting Metabolic Rate (Thermal Work Limit) as an Index of Thermal Stress. *Applied Occupational and Environmental Hygiene*, 17:3, 176–186. <https://doi.org/10.1080/104732202317070947>
