# GoSPL Configuration Examples

This directory contains example configuration files for GoSPL when used with DynEarthSol.

## Files

- `gospl_config_example.yml` - Basic example configuration showing typical parameters

## Configuration Parameters

The GoSPL configuration file uses YAML format with the following main sections:

### Domain
- `npx`, `npy`: Number of processors (usually 1 for DynEarthSol coupling)
- `minX`, `maxX`, `minY`, `maxY`: Domain bounds (km)
- `dx`, `dy`: Grid spacing (km)

### Time
- `start`, `end`: Start and end times (years)
- `dt`: Time step size (years)

### Surface Processes (`sp`)
- `hillslopeKa`: Hillslope diffusion coefficient (m²/yr)
- `streamPowerLaw`: Enable stream power law erosion
- `Kf`: Fluvial erosion coefficient
- `m`, `n`: Stream power law exponents
- `seaLevel`: Sea level (m)
- `marineKa`: Marine diffusion coefficient (m²/yr)

### Output
- `dir`: Output directory name
- `interval`: Output interval (years)

## Usage Notes

1. The domain bounds should match or encompass your DynEarthSol model domain
2. Time parameters are usually controlled by DynEarthSol, but need to be specified
3. Adjust erosion coefficients based on your geological setting
4. Output directory will be created relative to your DynEarthSol working directory
