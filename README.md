# Canadian Forest Fire Emission Prediction System #
Source code from Dr. Kerry Anderson of Canadian Forest Service, Natural Resources Canada, Government of Canada

Current Version of [CFFEPS User Manual](CFFEPS_20180607_v2.pdf)

The code posted here is the version of the CFFEPS source code that was used in the study described in the manuscript entitled **"The FireWork air quality forecast system with biomass burning emissions from the Canadian Forest Fire Emissions Prediction System"** that was submitted to the journal [Geoscientific Model Development](https://www.geoscientific-model-development.net) in 2019.

---

Once the files have been extracted, a file structure will be as follows:

    \Makefile
    \include    (contains all the include files required by CFFEPS)
    \obj        (contains the compiled objects)
    \src        (contains the source code written in C)
    \bin        (empty for compiled output)

In the \src directory will be the following file structure

    \src\CFFEPS     (contains the CFFEPS source code and relevant files)
    \src\diurnal    (contains diurnal adjustment calculations for the Fine Fuel Moisture Code FFFMC)
    \src\fbp2009    (contains the source code for the Canadaian Forest Fire Behaviour Prediction (FBP) System)
    \src\fwi84      (contains the source code for the Canadaian Forest Fire Weather Index (FWI) System)
    \src\thermo     (contains source code to do some rudimentary thermodynamic calculations)
    \src\utils      (conatins some utility functions for CFFEPS, such as reading and writing input files)

    
---

Internal revision tracking: this is cloned from arqi/cffeps repository (tag v2.03) revision 
*/cffeps/commit/6473bc036b96f8007f9d140bf684c52c17fdefde*
