# ocean_pole_tide.coef: station-code collisions give SMA and NOEMA the wrong coefficients (regression in e1957a846)

**Component**: difxcalc11: `data/ocean_pole_tide.coef` and the lookup in `src/dinit.f` (`dOPTLIN`)
**Introduced by**: commit `e1957a846` (2025-01-15, "update calc11 coeff files with
some VOGS/IVS stations from the most recent CDDIS blokq.c11.dat reference")
**Affected**: every build containing that commit (verified at v2.9.0-6).
Pre-2025 correlations (e.g. all EHT production) are unaffected, so current
DiFX silently disagrees with its own earlier production output.

## The problem

The CDDIS import added VGOS stations whose 2-char codes collide with
long-standing mm-VLBI codes, and dropped an alias line:

- `SW` now appears twice: `SMA` (Maunakea) at line 357 and the new
  `SVERT13V` (Svetloe, Russia, 60.5°N) at line 364. A station named `SW`
  resolves to **Svetloe**.
- The `NOEMA NN` alias line was **removed** (only `NOEMA NQ` remains) and
  `NN` is now carried by the new `NYALE13N` (Ny-Ålesund, Svalbard,
  78.9°N). A station named `NN` resolves to **Svalbard**; NOEMA's values
  are unreachable under that code.

DiFX/EHT correlator jobs name these telescopes by exactly these codes.

## Why it is silent

`dOPTLIN` (dinit.f) scans the whole file; every line matching a station by
8-char name **or** by 3-char code overwrites `OPTL6`; the last match wins:

```fortran
110   Read(iunit,1001,err=110,end=200) Vsite, C3, ...
      do I = 1, Numsit
        If (Vsite .eq. Dbsites(I) .or. C3 .eq. Dbsites(I)(1:3)) Then
          ! fills OPTL6(1..6,I); no exclusivity, no already-filled check
```

For `SW`, both lines 357 and 364 match; the later `SVERT13V` wins. The
wrong match also sets `Jsite`, so the "No ocean pole tide loading
coefficients" warning that would have exposed the problem never fires.

For contrast, `ocean_load.coef` kept consistent codes (`SVETLOE SV`,
`NOEMA NN` alias still present), and its one `NN` duplicate
(`NYALE13N` line 1666 vs `NOEMA` line 4141) gives the right answer only
because NOEMA happens to come later in the file. It is the same
last-match mechanism, correct by file order alone.

## Impact

Ocean pole tide loading is quasi-static over an observing track, so the
wrong coefficients produce a **constant sub-mm ITRF displacement** that
is indistinguishable from a station-position error. Measured against
EHT 2021 production correlator models (e21d15, DiFX 2.7.1, `.im` and
float64 SWIN output on the ALMA archive): SMA is off by ~0.3 mm and NOEMA
by ~0.1 mm, constant across sources and hours; with the two catalog lines
restored, both agree with production at the micrometre level. Any
current-DiFX correlation or model recomputation involving SMA or NOEMA
is affected.

## Suggested fix

Data (restores pre-`e1957a846` resolution):
1. Blank or change `SVERT13V`'s code `SW` (its full-name lookup still works).
2. Re-add the `NOEMA NN` alias line (same values as `NOEMA NQ`); blank
   `NYALE13N`'s `NN` or move it out of the mm-VLBI code space.
3. Audit both catalogs for the remaining duplicate codes
   (`ocean_pole_tide`: AT GB GO GT HD HH JB K2 KG MY N6 NT OO PT RO TA TI
   UD US WB YM; `ocean_load`: E1-E7 GB GL JB KT N1-N7 NN ON PT US W1-W6 WB)
   and resolve any that pair different physical sites.

Code (makes the class of bug impossible/loud): in `dOPTLIN`/`dOCNIN`,
give exact full-name matches precedence over code matches, take the first
match instead of the last, and warn when a station matches multiple lines
with differing values (also keeps `Jsite` honest so the missing-station
warning fires).
