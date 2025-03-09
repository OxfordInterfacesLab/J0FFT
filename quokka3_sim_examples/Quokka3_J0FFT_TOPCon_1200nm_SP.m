% example Quokka3 settingsfile
% simulates 2D carrier profile and 1D luminescence profile for a TOPCon cell with varying rear contact J0


% (c) 2025 Andreas Fell


SaveSpatial = 1;

Syntax = 'front and rear contact unit cell';

Solver.SolutionType = 'EL-PL';
Solver.ELPL.Type = 'PL';

Solver.Electrical.MetalModelType = 'constant-potential';
Solver.Electrical.SkinJ0nieff.Type = 'user-T';
Solver.Electrical.SkinJ0nieff.T = 25 degC;

Sweep.Enable = 1;
Sweep.NGroups = 1;
Sweep.GroupA(1).Parameter = 'FullRearSkin.Lumped.Electrical.ContactedRecombination.J0';
Sweep.GroupA(1).Values = lin[10,100,10];


UnitCellType = 'standard';
CellThickness = 150;
FrontMetalShape = 'H-pattern';
FrontContactShape = 'line';
FrontContactWidth = 3 um;
FrontFingerWidth = 3 um;
FrontContactPitch = 1.5 mm;
FrontFingerShadingFraction = 0;
FrontContact.OhmicResistivity = 1e-3;

RearContactShape = 'line';
RearMetalShape = 'H-pattern';
RearFingerWidth = 30 um;
RearContactWidth = 30 um;
RearContactPitch = 1.5 mm;
RearFingerShadingFraction = 1;
RearContact.OhmicResistivity = 1e-3;

FrontMetalPolarity = 'p-type';
RearMetalPolarity = 'n-type';

Thermal.T = 25 degC;

Bulk.BackgroundDoping.SettingType = 'dopingtype-resistivity';
Bulk.BackgroundDoping.DopingType = 'n-type';
Bulk.BackgroundDoping.Resistivity = 1.21;
Bulk.Electrical.Recombination.Type = 'intrinsic plus SRH';
Bulk.Electrical.Recombination.SRH(1).Type = 'tau-Et';
Bulk.Electrical.Recombination.SRH(1).Et_Ei = 0;
Bulk.Electrical.Recombination.SRH(1).taun = 0.5 ms;
Bulk.Electrical.Recombination.SRH(1).taup = 5 ms;
Bulk.Mesh.Quality = 'standard';

Optical.GenerationModelType = 'Text-Z';
Optical.FrontIllumination.Enable = 1;
Optical.FrontIllumination.Spectrum.Type = 'monochromatic';
Optical.FrontIllumination.Spectrum.MonoWavelength = 650 nm;
Optical.FrontIllumination.Spectrum.MonoFlux = 3.275e17;
Optical.FrontIllumination.Scale = 1;
Optical.RearIllumination.Enable = 0;
Optical.TextZ.FrontText.Type = 'Text';
Optical.TextZ.FrontText.Text.Type = 'const';
Optical.TextZ.FrontText.Text.Const = 0.93;
Optical.TextZ.FrontText.FacetAngle = 54;
Optical.TextZ.FrontZ.Type = 'Basore';
Optical.TextZ.Basore.Front.Rint = 0.92;
Optical.TextZ.Basore.Front.Diffuse = 1;
Optical.TextZ.Basore.Rear.Rint = 0.92;
Optical.TextZ.Basore.Rear.Appp = 0.03;
Optical.TextZ.Basore.Rear.Diffuse = 1;

% defines spectral transmission of the luminescence detection optics
% sharp cut-off at 1200nm represents a perfect 1200nm short-pass filter
Optical.Luminescence.OpticsTransmission.Type = 'XY';
Optical.Luminescence.OpticsTransmission.X = [100 1200 1201 2000];
Optical.Luminescence.OpticsTransmission.Y = [1 1 0 0];


% full area front skin properties (p+ emitter)
FullFrontSkin.ElectricalModelType = 'lumped';
FullFrontSkin.Lumped.Electrical.RsheetEnable = 1;
FullFrontSkin.Lumped.Electrical.Rsheet = 150;
FullFrontSkin.Lumped.Electrical.NonContactedRecombination.ModelType = 'J0';
FullFrontSkin.Lumped.Electrical.NonContactedRecombination.J0 = 25 fA/cm2;
FullFrontSkin.Lumped.Electrical.NonContactedRecombination.J02 = 1 nA/cm2;
FullFrontSkin.Lumped.Electrical.ContactedRecombination.ModelType = 'J0';
FullFrontSkin.Lumped.Electrical.ContactedRecombination.J0 = 25 fA/cm2;
FullFrontSkin.Lumped.Electrical.ContactedRecombination.J02 = 1 nA/cm2;


% full rear skin properties (TOPCon)
FullRearSkin.ElectricalModelType = 'lumped';
FullRearSkin.Lumped.Electrical.RsheetEnable = 1;
FullRearSkin.Lumped.Electrical.Rsheet = 50;
FullRearSkin.Lumped.Electrical.NonContactedRecombination.ModelType = 'J0';
FullRearSkin.Lumped.Electrical.NonContactedRecombination.J0 = 10 fA/cm2;
FullRearSkin.Lumped.Electrical.ContactedRecombination.ModelType = 'J0';
FullRearSkin.Lumped.Electrical.ContactedRecombination.J0 = 50 fA/cm2;



