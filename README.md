# Code and data for "Layered national and provincial carbon policies in Sichuan and Chongqing: a small local temperature signal and little near-term heatwave relief"

- `Research_package/` contains the code, inputs, ensemble results and tables; see `Research_package/README.md` for the contents and how to run the analysis.
- The two daily ERA5-Land temperature arrays (`inputs/weather/tmax_float64.npy` and `tmin_float64.npy`, 545 MB each) exceed GitHub's file-size limit and are stored in `parts/` as numbered parts of at most 10 MB. Run `python parts/join_parts.py`, then extract `Research_package_tmax.zip` and `Research_package_tmin.zip` into the repository root; this places the arrays in `Research_package/inputs/weather/`. They are needed only for the `weather` stage.
- `parts/` also contains the parts of `Research_package_main.zip`, an archive of the same files as `Research_package/`.
