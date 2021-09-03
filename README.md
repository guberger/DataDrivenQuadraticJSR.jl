# DataDrivenQuadraticJSR.jl

## Workflow

### To run the GUI
```julia
(@v1.6) pkg> activate ./examples
```
Run the `main.jl` example by executing the corresponding file, or running
```julia
julia> include("./examples/main.jl")
```
![GUI](https://github.com/guberger/DataDrivenQuadraticJSR.jl/blob/main/gui_picture.PNG)

You can fine-tune (even go out of bounds) the value of the sliders using the following commands in the REPL:
```julia
julia> my_gui.slider_γ.set_val(0.889319) # choose your val
julia> my_gui.slider_tr.set_val(2.319949) # choose your val
julia> my_gui.slider_ηm.set_val(10.0) # choose your val
julia> my_gui.slider_ηa.set_val(1e-5) # choose your val
```
