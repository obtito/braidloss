/*
 * Q1_SolidWire_generated.java
 */

import com.comsol.model.*;
import com.comsol.model.util.*;

/** Model exported on Sep 25 2026, 19:36 by COMSOL 6.2.0.290. */
public class Q1_SolidWire_generated {

  public static Model run() {
    Model model = ModelUtil.create("Model");

    model.modelPath("C:\\Users\\zq257\\Documents\\Codex\\2026-09-25\\jin\\outputs\\Q1_COMSOL");

    model.label("Q1 - Solid copper wire - 200 kHz - 20 A RMS");

    model.param().set("a", "1[mm]", "Copper radius");
    model.param().set("L", "1[m]", "Evaluated wire length");
    model.param().set("sigmaCu", "5.8e7[S/m]", "Copper conductivity at 20 degC");
    model.param().set("muCu", "4*pi*1e-7[H/m]", "Permeability (mu_r=1)");
    model.param().set("Irms", "20[A]", "RMS transport current");
    model.param().set("Ipk", "sqrt(2)*Irms", "Peak phasor current for frequency domain");
    model.param().set("f0", "200[kHz]", "Operating frequency");
    model.param().set("delta0", "sqrt(1/(pi*f0*muCu*sigmaCu))", "Theoretical skin depth");
    model.param().set("Rdc", "L/(sigmaCu*pi*a^2)", "DC resistance");
    model.param().set("Nr", "16", "Radial mapped elements");
    model.param().set("Nz", "4", "Axial mapped elements (z invariant solution)");

    model.component().create("comp1", true);

    model.component("comp1").geom().create("geom1", 2);
    model.component("comp1").geom("geom1").axisymmetric(true);
    model.component("comp1").geom("geom1").create("r1", "Rectangle");
    model.component("comp1").geom("geom1").feature("r1").set("size", new String[]{"a", "L"});
    model.component("comp1").geom("geom1").run();

    model.component("comp1").material().create("mat1", "Common");
    model.component("comp1").material("mat1").label("Oxygen-free copper, 20 degC");
    model.component("comp1").material("mat1").propertyGroup("def")
         .set("electricconductivity", new String[]{"sigmaCu"});
    model.component("comp1").material("mat1").propertyGroup("def")
         .set("relpermeability", new String[]{"muCu/mu0_const"});
    model.component("comp1").material("mat1").propertyGroup("def").set("relpermittivity", new String[]{"1"});

    model.component("comp1").physics().create("mfh", "MagneticFieldFormulation", "geom1");
    model.component("comp1").physics("mfh").label("H-phi formulation: axial transport current");
    model.component("comp1").physics("mfh").prop("components").set("components", "inplane");
    model.component("comp1").physics("mfh").prop("DivergenceConstraint").set("DivergenceConstraint", false);
    model.component("comp1").physics("mfh").create("mfb1", "MagneticFieldBoundary", 1);
    model.component("comp1").physics("mfh").feature("mfb1").selection().set(4);
    model.component("comp1").physics("mfh").feature("mfb1").label("Surface Hphi = Ipeak / (2 pi a)");
    model.component("comp1").physics("mfh").feature("mfb1").set("H0", new String[]{"0", "Ipk/(2*pi*a)", "0"});
    model.component("comp1").physics("mfh").create("axisH", "PointwiseConstraint", 1);
    model.component("comp1").physics("mfh").feature("axisH").selection().set(1);
    model.component("comp1").physics("mfh").feature("axisH").label("Axis regularity: Hphi = 0 at r = 0");
    model.component("comp1").physics("mfh").feature("axisH").set("constraintExpression", "Hphi");
    model.component("comp1").physics("mfh").feature("fl1").featureInfo("info")
         .set("mfh.curlHz", new String[]{"Hphir+if(abs(r)<1e-12[m],Hphir,Hphi/r)"});

    model.component("comp1").cpl().create("intV", "Integration");
    model.component("comp1").cpl("intV").selection().all();

    model.component("comp1").variable().create("var1");
    model.component("comp1").variable("var1").set("Jrms", "mfh.normJ/sqrt(2)", "RMS magnitude of current density");
    model.component("comp1").variable("var1")
         .set("qloss", "mfh.normJ^2/(2*sigmaCu)", "Time-average Joule loss density");

    model.component("comp1").cpl("intV").label("Meridian area integral; use explicit 2*pi*r Jacobian");

    model.component("comp1").variable("var1").set("Pcu", "intV(2*pi*r*qloss)", "Copper loss over 1 m");
    model.component("comp1").variable("var1").set("Rac", "Pcu/Irms^2", "AC resistance from Joule losses");
    model.component("comp1").variable("var1").set("ratio", "Rac/Rdc", "AC/DC resistance ratio");
    model.component("comp1").variable("var1")
         .set("Icheck", "intV(2*pi*r*mfh.Jz)/L/sqrt(2)", "Complex RMS current, axial average");

    model.component("comp1").mesh().create("mesh1");
    model.component("comp1").mesh("mesh1").create("map1", "Map");
    model.component("comp1").mesh("mesh1").feature("map1").create("disr", "Distribution");
    model.component("comp1").mesh("mesh1").feature("map1").feature("disr").selection().set(2, 3);
    model.component("comp1").mesh("mesh1").feature("map1").feature("disr").set("numelem", "Nr");
    model.component("comp1").mesh("mesh1").feature("map1").create("disz", "Distribution");
    model.component("comp1").mesh("mesh1").feature("map1").feature("disz").selection().set(1, 4);
    model.component("comp1").mesh("mesh1").feature("map1").feature("disz").set("numelem", "Nz");

    model.study().create("std1");
    model.study("std1").label("200 kHz frequency-domain solution");
    model.study("std1").create("freq", "Frequency");
    model.study("std1").feature("freq").set("plist", "f0");

    model.result().numerical().create("gev", "EvalGlobal");
    model.result().numerical("gev")
         .set("expr", new String[]{"freq", "delta0", "Rdc", "comp1.Rac", "comp1.ratio", "comp1.Pcu", "real(comp1.Icheck)", "imag(comp1.Icheck)"});
    model.result().numerical("gev").set("unit", new String[]{"Hz", "m", "ohm", "ohm", "1", "W", "A", "A"});
    model.result().table().create("tbl1", "Table");
    model.result().table("tbl1").label("Electrical results (RMS convention)");
    model.result().numerical("gev").set("table", "tbl1");

    model.param().set("Nr", "8");

    model.component("comp1").mesh("mesh1").run();

    model.study("std1").run();

    model.result().numerical("gev").set("data", "dset1");
    model.result().numerical("gev").setResult();

    model.param().set("Nr", "16");

    model.component("comp1").mesh("mesh1").run();

    model.study("std1").run();

    model.result().numerical("gev").set("data", "dset1");
    model.result().numerical("gev").setResult();

    model.param().set("Nr", "32");

    model.component("comp1").mesh("mesh1").run();

    model.study("std1").run();

    model.result().numerical("gev").set("data", "dset1");
    model.result().numerical("gev").setResult();

    model.param().set("Nr", "64");

    model.component("comp1").mesh("mesh1").run();

    model.study("std1").run();

    model.result().numerical("gev").set("data", "dset1");
    model.result().numerical("gev").setResult();

    model.param().set("Nr", "128");

    model.component("comp1").mesh("mesh1").run();

    model.study("std1").run();

    model.result().numerical("gev").set("data", "dset1");
    model.result().numerical("gev").setResult();
    model.result().numerical().create("rad", "Interp");
    model.result().numerical("rad").set("data", "dset1");
    model.result().numerical("rad")
         .set("expr", new String[]{"r", "z", "Jrms", "real(mfh.Jz)/sqrt(2)", "imag(mfh.Jz)/sqrt(2)", "qloss", "abs(mfh.Jr)/sqrt(2)"});

    model.label("Q1_SolidWire_200kHz.mph");

    model.result().dataset().create("radline", "CutLine2D");
    model.result().dataset("radline").label("Radius at wire midspan z = 0.5 m");
    model.result().dataset("radline").set("data", "dset1");
    model.result().dataset("radline").set("genpoints", new double[][]{{0, 0.5}, {0.001, 0.5}});
    model.result().create("pgRad", "PlotGroup1D");
    model.result("pgRad").label("Radial current density - RMS");
    model.result("pgRad").set("data", "radline");
    model.result("pgRad").create("line1", "LineGraph");
    model.result("pgRad").feature("line1").set("expr", "Jrms");
    model.result("pgRad").feature("line1").set("unit", "A/mm^2");
    model.result("pgRad").feature("line1").set("xdata", "expr");
    model.result("pgRad").feature("line1").set("xdataexpr", "r");
    model.result("pgRad").feature("line1").set("xdataunit", "mm");
    model.result("pgRad").feature("line1").set("resolution", "extrafine");
    model.result().dataset().create("revWire", "Revolve2D");
    model.result().dataset("revWire").set("data", "dset1");
    model.result().dataset("revWire").set("revangle", 360);
    model.result().dataset("revWire").set("layermethod", "custom");
    model.result().dataset("revWire").set("revlayers", 160);
    model.result().dataset("revWire").label("Full cylindrical wire, D=2 mm, L=1 m");
    model.result().dataset().create("cross", "CutPlane");
    model.result().dataset("cross").set("data", "revWire");
    model.result().dataset("cross").set("quickplane", "xy");
    model.result().dataset("cross").set("quickz", "L/2");
    model.result().dataset("cross").label("Wire transverse section at z = L/2");
    model.result().create("pgCross", "PlotGroup2D");
    model.result("pgCross").label("Cross section - current density RMS (A/mm^2)");
    model.result("pgCross").set("data", "cross");
    model.result("pgCross").create("surf1", "Surface");
    model.result("pgCross").feature("surf1").set("expr", "Jrms");
    model.result("pgCross").feature("surf1").set("unit", "A/mm^2");
    model.result("pgCross").feature("surf1").set("colortable", "ThermalLight");
    model.result("pgCross").feature("surf1").set("resolution", "fine");
    model.result().create("pgWire", "PlotGroup3D");
    model.result("pgWire").label("Full 1 m copper wire (true aspect ratio)");
    model.result("pgWire").set("data", "revWire");
    model.result("pgWire").create("surf1", "Surface");
    model.result("pgWire").feature("surf1").set("expr", "Jrms");
    model.result("pgWire").feature("surf1").set("unit", "A/mm^2");
    model.result("pgWire").feature("surf1").set("colortable", "ThermalLight");
    model.result("pgCross").run();
    model.result("pgRad").run();
    model.result().table().create("tblMesh", "Table");
    model.result().table("tblMesh").label("Mesh convergence: 8 / 16 / 32 / 64 / 128 radial elements");
    model.result().table("tblMesh")
         .importData("C:/Users/zq257/Documents/Codex/2026-09-25/jin/outputs/Q1_COMSOL/mesh_convergence.csv");
    model.result().numerical("gev").setResult();

    model.label("Q1_SolidWire_200kHz.mph");

    model.result("pgCross").set("titletype", "manual");
    model.result("pgCross").set("title", "200 kHz | Current density RMS (A/mm^2) | z = 0.5 m");
    model.result().export().create("imgCross", "pgCross", "Image");
    model.result().export("imgCross")
         .set("pngfilename", "C:/Users/zq257/Documents/Codex/2026-09-25/jin/outputs/Q1_COMSOL/comsol_cross_section.png");
    model.result().export("imgCross").set("size", "manualweb");
    model.result().export("imgCross").set("unit", "px");
    model.result().export("imgCross").set("width", 1500);
    model.result().export("imgCross").set("height", 1100);
    model.result().export("imgCross").set("options2d", "on");
    model.result().export("imgCross").set("legend2d", "on");
    model.result().export("imgCross").set("axes2d", "on");
    model.result().export("imgCross").set("title2d", "on");
    model.result().export("imgCross").set("zoomextents", true);
    model.result().export("imgCross").run();
    model.result("pgCross").run();

    return model;
  }

  public static void main(String[] args) {
    run();
  }

}
