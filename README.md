# Normal Depth in Trapezoidal Open Channels by the Newton–Raphson Method (MATLAB)

CEE 215: Numerical Analysis and Computer Programming, Section 1, North South University (Fall 2026).

Manning's equation, `Q = (k/n) A R^(2/3) S0^(1/2)`, cannot be solved for the flow depth directly. This program
rewrites it as `f(y) = 0` and finds the **normal depth** with the Newton–Raphson method, using an
analytical derivative `f'(y)`. The same solver also gives the **critical depth**. The program then reports
the velocity, Froude number, flow regime and channel slope type, and checks every answer against MATLAB's `fzero`.

## Files

| File | Purpose |
|---|---|
| `matlab/NR_Channel_Main.m` | **Main program, run this.** Menu: Example 1 (SI canal), Example 2 (Chow 1959 Ex. 6.2), or your own data |
| `matlab/newton_raphson.m` | General Newton–Raphson solver: root, iterations, history table, status flag |
| `matlab/read_input.m` | Keyboard input with a default value and validity checks |
| `matlab/NR_Limitations_Demo.m` | Eight worked examples of the limitations (zero/negative guess, S0 = 0, oscillation, divergence, multiple root, circular culvert) |
| `matlab/show_case.m` | Prints iteration tables for the demo |
| `matlab/results/` | Figures and `NR_Channel_Results.txt` written by the scripts |
| `docs/Proposal_CEE215_NewtonRaphson.docx / .pdf` | Project proposal (4 pages, with flowchart) |
| `docs/Final_Report_CEE215_NewtonRaphson.docx / .pdf` | Final report (theory, flowchart, code, results, limitations, contribution chart) |
| `docs/figures/` | Flowchart and figures used in the documents |

## How to run

1. Open MATLAB and set the `matlab` folder as the Current Folder.
2. Type `NR_Channel_Main` and press Enter. Choose 1, 2 or 3. Pressing Enter at any prompt accepts the default value shown in brackets.
3. Type `NR_Limitations_Demo` to step through the limitation examples. Press Enter between cases.

The code also runs in GNU Octave.

## Key results

| Case | Normal depth yn | Iterations | Check |
|---|---|---|---|
| Example 1: b = 6 m, z = 1.5, n = 0.025, S0 = 0.0004, Q = 25 m³/s | 2.403895 m | 4 | `fzero` agrees to 1e-13 % |
| Example 2: Chow (1959) Ex. 6.2, b = 20 ft, z = 2, n = 0.025, S0 = 0.0016, Q = 400 cfs | 3.360968 ft | 4 | Published answer 3.36 ft |
