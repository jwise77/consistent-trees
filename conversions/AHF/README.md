AHF to CTrees Coversion
=======================

The script was written by [Yao-Yuan Mao](https://yymao.github.io) for the [Sussing Merger Trees](https://arxiv.org/find/astro-ph/1/ti:+EXACT+Sussing_Merger_Trees/0/1/0/all/0/1) comparison project.

## Usage

The `convert_input.py` script simply rearranges the columns. The complication comes in as Rockstar prints out `DescID` while AHF does not (at least at the time when I wrote the script). Hence, `get_descendants.c` is needed to generate `DescID`, which will later be read in by `convert_input.py`.

1. Run `make` to compile `get_descendants.c`

2. Run `./get_descendants <particle_list_snap1> <particle_list_snap2> <particle_list_snap3> ...`, where each of the particle list is from one epoch (snapshot). The order of the particle lists should be from high to low redshifts. See the next section for the assumed format of the particle list.

3. Run `python convert_input.py <DIR_AHF_HALOS> <DIR_DESC_BIN>`, where `DIR_AHF_HALOS` and `DIR_DESC_BIN` are the directories that contain `*.AHF_halos` and `desc_*.bin` files respectively.

## Format of the particle lists

```
<N = Number of halos>
<ID if Halo 1> <M1 = Number of particles of Halo 1>
<ID of Particle 1 of Halo 1>
<ID of Particle 2 of Halo 1>
...
<ID of Particle M1 of Halo 1>
<ID if Halo 2> <M2 = Number of particles of Halo 2>
<ID of Particle 1 of Halo 2>
<ID of Particle 2 of Halo 2>
...
<ID of Particle M2 of Halo 2>
...
<ID if Halo N> <MN = Number of particles of Halo N>
<ID of Particle 1 of Halo N>
<ID of Particle 2 of Halo N>
...
<ID of Particle MN of Halo N>
```
