# multimodal-annotation-distance
An [ELAN](https://archive.mpi.nl/tla/elan) Python extension for determining distances between multimodal annotations.

## Dependencies
- [pimpy](https://github.com/dopefishh/pympi)
- [pandas](https://pandas.pydata.org/)
- [tqdm](https://github.com/tqdm/tqdm)

## Installation
1. Clone this repository on your local computer.
2. In the command line, go to the repository containing the `multimodal_annotation_distance.py` file.
3. Create and activate a virtual environment of your preference (e.g., [Conda](https://anaconda.org/anaconda/conda) or [virtualenv](https://virtualenv.pypa.io/en/latest/)).
4. To install the required dependencies, you can use the provided `requirements.txt` file by typing `pip install -r requirements.txt` in the command line.

## Usage
1. Activate your Python virtual environment.
2. If  you haven't done so before, you need to fill in the required values inside the provided `config.ini` file. Take care of not erasing any of them from the file!
3. To run the program, inside the directory where the `multimodal_annotation_distance.py` file is located, type the following in the command line: `python multimodal_annotation_distance.py --config-file-path path/to/config.ini`
