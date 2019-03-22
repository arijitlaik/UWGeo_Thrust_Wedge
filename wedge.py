# %% markdown
# # Thust Wedge Tutorial
#
# Model Setup: Arijit Laik(a.laik@vu.nl)
#
# Revision: Romain Beucher (rbeucher@unimelb.edu.au)
#
# This tutorial is mostly based off Ruh et al. (2013) and other similar stuides, which examine the development thrust sheet and or accretionary wegdes with visco-plastic rheologies .In this 2D model, the boundary conditions, initial conditions,geometry material properties are based on Ruh et al. and other such similar studies.(see Ruh et al. 2013 and references therein: J.B. Ruh, T. Gerya nd J.-P. Burg (2013), G3, v.14(4), p. 1131-1155)
# %% markdown
# ![Tutorial10](./images/Tutorial_10.gif)
# %%
import underworld.function as fn
import UWGeodynamics as GEO
import numpy as np
import glucifer

# %% markdown
# ### Model Scaling
# %%
u = GEO.UnitRegistry

velocity = 1 * u.centimeter / u.year
model_length = 100.0 * u.kilometer
bodyforce = 2700.0 * u.kilogram / u.metre ** 3 * 9.81 * u.meter / u.second ** 2

KL = model_length
Kt = KL / velocity
KM = bodyforce * KL ** 2 * Kt ** 2

GEO.scaling["[length]"] = KL
GEO.scaling["[time]"] = Kt
GEO.scaling["[mass]"] = KM

Model = GEO.Model(
    elementRes=(512, 256),
    minCoord=(0.0 * u.kilometer, -7 * u.kilometer),
    maxCoord=(128.0 * u.kilometer, 9.0 * u.kilometer),
    gravity=(0.0, -9.81 * u.meter / u.second ** 2),
)
# %%
Model.outputDir = "outputs_tutorial10"
# %% markdown
# ### Material Setup
#
# We will start by defining the air layer, the rigid base and the frictional layer.
# %%
air_shape = GEO.shapes.Layer(top=Model.top, bottom=0.0 * u.kilometer)
fricLayerShape = GEO.shapes.Layer(
    top=Model.bottom + 1.0 * u.kilometer, bottom=Model.bottom + 0.5 * u.kilometer
)
rigidBaseShape = GEO.shapes.Layer(
    top=Model.bottom + 0.5 * u.kilometer, bottom=Model.bottom
)

air = Model.add_material(name="Air", shape=air_shape)
frictionalBasal = Model.add_material(name="Frictional", shape=fricLayerShape)
rigidBase = Model.add_material(name="Frictional", shape=rigidBaseShape)
sediment = Model.add_material(name="Sediment")
# %% markdown
# ### Pile of sedimentary layers
# %% markdown
# We define a 6km thick pile of sedimentary layers, with each layer being 0.5 km thick, for a total of 12 layers
# %%
top_pile = 0.0
bottom_pile = -6.0 * u.kilometer

NLayers = 12
layer_thickness = (top_pile - bottom_pile) / NLayers

plastic_pile = []

layer_above = air_shape

for index in range(NLayers):
    shape = GEO.shapes.Layer(
        top=layer_above.bottom, bottom=layer_above.bottom - layer_thickness
    )
    material = Model.add_material(name="Plastic {0}".format(index), shape=shape)
    plastic_pile.append(material)
    layer_above = shape
# %%
colours = "tan #3A3A3A #00a8a8 " + 6 * "#425e6a salmon "

FigMat = glucifer.Figure(figsize=(1200, 250), quality=3)
FigMat.Points(
    Model.swarm,
    Model.materialField,
    discrete=True,
    colours=colours,
    fn_mask=Model.materialField > air.index,
)
FigMat.objects[0].colourBar["binlabels"] = True
FigMat.show()
# %% markdown
# ## Model (Global) properties
# %%
Model.density = 2700 * u.kilogram / u.metre ** 3
Model.viscosity = 1e23 * u.pascal * u.second
Model.maxViscosity = 1e23 * u.pascal * u.second
Model.minViscosity = 5e19 * u.pascal * u.second
# %% markdown
# ## Viscosity
# %%
air.viscosity = 1e19 * u.pascal * u.second
air.minViscosity = 1e19 * u.pascal * u.second
air.density = 1.0 * u.kilogram / u.metre ** 3

# Note that this is not necessary as this does not differ from the
# Model property.
for material in plastic_pile:
    material.density = 2700 * u.kilogram / u.metre ** 3
    material.viscosity = 1e23 * u.pascal * u.second

frictionalBasal.viscosity = 1e23 * u.pascal * u.second
rigidBase.viscosity = 1e23 * u.pascal * u.second
# %% markdown
# ## Plasticity
# %%
plastic_Law = GEO.DruckerPrager(
    cohesion=20.0 * u.megapascal,
    cohesionAfterSoftening=4.0 * u.megapascal,
    frictionCoefficient=np.tan(np.radians(25.0)),
    frictionAfterSoftening=np.tan(np.radians(20.0)),
    epsilon1=0.01,
    epsilon2=0.06,
)
sediment.plasticity = GEO.DruckerPrager(
    cohesion=10.0 * u.megapascal,
    cohesionAfterSoftening=4.0 * u.megapascal,
    frictionCoefficient=np.tan(np.radians(20.0)),
    frictionAfterSoftening=np.tan(np.radians(15.0)),
    epsilon1=0.01,
    epsilon2=0.06,
)
for material in plastic_pile:
    material.plasticity = plastic_Law

frictionalBasal.plasticity = GEO.DruckerPrager(
    cohesion=0.1 * u.megapascal,
    frictionCoefficient=np.tan(np.radians(12.0)),
    frictionAfterSoftening=np.tan(np.radians(6.0)),
    epsilon1=0.01,
    epsilon2=0.06,
)
# %% markdown
# ## Velocity Boundary Conditions
# The rigid base has the same velocity(=1cm/year,or velocity of the left wall), the linear gradational velocity is in the low frictional layer, this prevents the rigid bottom sheet from bending at the left edge, moreover as the rigid bottom sheet is analogous to a Mylar(PET) conveyor belt / sheet in sandbox models(Konstantinovskaya and Malavieille, 2011[https://doi.org/10.1016/j.tecto.2011.01.020], or Bose et al 2014[https://doi.org/10.1016/j.jsg.2014.07.004])
# The boundary conditions simulate the mechanics of convergent plate boundaries where a rigid backstop scrapes upper levels of the crust off a rigid moving "plate". This is comparable to setups of analog models, where a rough sheet lying below sand layers is pulled out below a fixed and rigid backstop .This setup also matches the boundary conditions used in the analytical critical wedge theory. New swarms coming in through the left side allow for simulation of long term deformation.
#
# ![setupandBCs](./images/Tutorial_10_bcs.png)

# %%

tapeL = frictionalBasal
flthick = GEO.nd(tapeL.top - tapeL.bottom)

conditions = [
    (Model.y <= GEO.nd(rigidBase.top), GEO.nd(-velocity)),
    (
        Model.y < GEO.nd(tapeL.top),
        GEO.nd(-velocity) * (flthick - (Model.y - GEO.nd(tapeL.bottom))) / flthick,
    ),
    (True, GEO.nd(0.0 * u.centimeter / u.year)),
]

fn_condition = fn.branching.conditional(conditions)

Model.set_velocityBCs(
    left=[fn_condition, 0.0],
    right=[-velocity, None],
    top=[None, None],
    bottom=[-velocity, 0.0],
)
# %%
Fig = glucifer.Figure(figsize=(1200, 250), quality=3)
Fig.Surface(Model.mesh, fn.math.dot(Model.velocityField, Model.velocityField))
Fig.show()
# %%
Model.solver.set_inner_method("mumps")
Model.solver.set_penalty(1e6)
GEO.rcParams["nonlinear.tolerance"] = 1e-3
GEO.rcParams["initial.nonlinear.tolerance"] = 1e-5
# %%
Model.init_model()
# %%
Fig = glucifer.Figure(figsize=(1200, 250), quality=3)
Fig.Points(
    Model.swarm,
    GEO.Dimensionalize(Model.viscosityField, u.pascal * u.second),
    logScale=True,
)
Fig.show()
# %%
Model.surfaceProcesses = GEO.surfaceProcesses.Badlands(
    airIndex=[air.index],
    sedimentIndex=sediment.index,
    XML="ressources/badlandsT10.xml",
    resolution=2.0 * u.kilometer,
    checkpoint_interval=0.01 * u.megayears,
)
# %%


def post_solve_hook():
    global FigMat
    if Model.step % 10 == 0:
        FigMat.save("Material-{0}.png".format(Model.step))


Model.postSolveHook = post_solve_hook
# %%
Model.run_for(0.2 * u.megayears, checkpoint_interval=0.1 * u.megayear, restartStep=None)
# %%

# %%
air.index
# %%
