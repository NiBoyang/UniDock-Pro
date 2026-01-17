# UniDock‑Pro

**UniDock‑Pro** is an all‑in‑one **GPU‑accelerated** virtual screening platform. It provides **classical docking**, **ligand similarity searching**, and **hybrid docking** in a single binary.

> **Project lineage.** UniDock‑Pro is a fork of **[Uni‑Dock](https://github.com/dptech-corp/Uni-Dock/)**, which is a GPU-accelerated version of **[AutoDock-Vina](https://github.com/ccsb-scripps/AutoDock-Vina)**. Please also see the Uni‑Dock publication: **DOI: 10.1021/acs.jctc.2c01145**.

## Citation & Preprint

If you use UniDock‑Pro in academic work, please cite the Uni‑Dock paper above and our preprint:

* **UniDock‑Pro preprint:** [https://chemrxiv.org/engage/chemrxiv/article-details/689e7544728bf9025e86ce58](https://chemrxiv.org/engage/chemrxiv/article-details/689e7544728bf9025e86ce58)

---

## Building from source

1. **Install dependencies**

   * **CUDA Toolkit ≥ 11.8** — follow NVIDIA’s [installation guide](https://docs.nvidia.com/cuda/cuda-installation-guide-linux/index.html).
   * **CMake ≥ 3.16**
   * **C++ compiler** (must be [compatible with NVCC](https://docs.nvidia.com/cuda/cuda-installation-guide-linux/index.html#host-compiler-support-policy); `g++` works in most cases)
   * **Boost ≥ 1.72**

     ```bash
     sudo apt install \
       libboost-system-dev libboost-thread-dev libboost-serialization-dev \
       libboost-filesystem-dev libboost-program-options-dev libboost-timer-dev
     ```

     If the above does not satisfy version requirements, install from [Boost source](https://www.boost.org/users/download/) or via Anaconda:
     `conda install -c anaconda libboost`.

2. **Clone the repository**

   ```bash
   git clone https://github.com/NiBoyang/UniDock-Pro.git
   ```

3. **Build with CMake**

   ```bash
   cd UniDock-Pro
   cmake -B build
   cmake --build build -j4
   ```

The main binary `udp` will be generated under `build/`.

---

## Usage

Run `./udp --help` for the full set of options. Minimal, reproducible examples are shown below.

> **Required option — `--search_mode`.**
> You **must** explicitly set `--search_mode` to one of `fast`, `balance`, or `detail`. Leaving it unset is currently **unsupported** and may lead to **incorrect or unstable behavior**.

### Input expectations

* **Receptor:** `.pdbqt`
* **Ligands:** Two input methods are supported:
  * A text file (`--ligand_index`) listing `.pdbqt` paths, one per line
  * A directory (`--ligand_directory`) containing `.pdbqt` files
* **Reference ligand(s) (for similarity searching or hybrid docking):**
  * Single reference: `--reference_ligand ref.pdbqt`
  * Multiple references (ensemble): `--reference_ligand ref1.pdbqt ref2.pdbqt ref3.pdbqt` (space-separated)
* **Search box:** `--center_{x,y,z}` and `--size_{x,y,z}` in Å
* **Output directory (optional):** `--dir <path>`

### Classical docking

Using a ligand index file:

```bash
./udp \
  --receptor rec.pdbqt \
  --ligand_index ligand_index.txt \
  --center_x 0 --center_y 0 --center_z 0 \
  --size_x 20 --size_y 20 --size_z 20 \
  --search_mode balance \
  --dir ./results
```

Using a ligand directory:

```bash
./udp \
  --receptor rec.pdbqt \
  --ligand_directory ./ligands/ \
  --center_x 0 --center_y 0 --center_z 0 \
  --size_x 20 --size_y 20 --size_z 20 \
  --search_mode balance \
  --dir ./results
```

### Ligand similarity searching

With a single reference ligand:

```bash
./udp \
  --reference_ligand ref.pdbqt \
  --ligand_index ligand_index.txt \
  --center_x 0 --center_y 0 --center_z 0 \
  --size_x 20 --size_y 20 --size_z 20 \
  --search_mode fast \
  --dir ./results
```

With multiple reference ligands (ensemble):

```bash
./udp \
  --reference_ligand ref1.pdbqt ref2.pdbqt ref3.pdbqt \
  --ligand_index ligand_index.txt \
  --center_x 0 --center_y 0 --center_z 0 \
  --size_x 20 --size_y 20 --size_z 20 \
  --search_mode fast \
  --dir ./results
```

### Hybrid docking

> **Important:** In **hybrid mode**, the `--reference_ligand` **must be provided in its co‑crystallized pose with the specified receptor**.

With a single reference ligand:

```bash
./udp \
  --receptor rec.pdbqt \
  --reference_ligand ref.pdbqt \
  --ligand_index ligand_index.txt \
  --center_x 0 --center_y 0 --center_z 0 \
  --size_x 20 --size_y 20 --size_z 20 \
  --search_mode detail \
  --dir ./results
```

With multiple reference ligands (ensemble):

```bash
./udp \
  --receptor rec.pdbqt \
  --reference_ligand ref1.pdbqt ref2.pdbqt ref3.pdbqt \
  --ligand_index ligand_index.txt \
  --center_x 0 --center_y 0 --center_z 0 \
  --size_x 20 --size_y 20 --size_z 20 \
  --search_mode detail \
  --dir ./results
```

### Worked example (from the `example/` folder)

From the repository’s `example/` directory:

```bash
  ../build/udp \
  --receptor ./receptor/rec.pdbqt \
  --reference_ligand ./ref_lig/xtal_lig.pdbqt \
  --ligand_index ligand_index.txt \
  --center_x 32.790 --center_y 38.342 --center_z 58.486 \
  --size_x 28 --size_y 28 --size_z 28 \
  --search_mode balance \
  --dir ./results
```

---

## Acknowledgements

UniDock‑Pro builds on the ideas and engineering of **AutoDock Vina** and the GPU‑accelerated **Uni‑Dock** project. We thank their authors and contributors for foundational work that enabled this fork.