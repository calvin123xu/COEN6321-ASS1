﻿/* Qiaoyu Xu 40292010 */
// pls put the input file in the root folder which is the same folder where mian.cpp is in
// the output file will also generated in that folder

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <iomanip>    // For formatting output
#include <cstdlib>    // For rand() and srand()
#include <ctime>      // For time()
#include <unordered_map>
#include <algorithm>  // For std::min_element and std::unique // For std::sort
#include <utility>  // For std::pair
#include <unordered_set>
#include <functional>
#include <random>
#include <cmath>        // For exp()
#include <limits>  // for std::numeric_limits

using namespace std;

// a struct to represent each individual with its traids
struct Individual {
	vector<vector<string>> puzzle;
	int fitness;
	int age;// age is in the unit of generation, age = 2 mean it has been there for two generation
};

// for test purpose
void printVector(const std::vector<std::vector<std::string>>& vec) {
	for (const auto& row : vec) {         // Iterate over each row
		for (const auto& element : row) {  // Iterate over each element in the row
			std::cout << element << " ";    // Print the element followed by a space
		}
		std::cout << std::endl;             // Print a newline after each row
	}
}


// this method simply reads the file Ass1Input in the root directory
vector<vector<string>> readFile() {
	ifstream inputFile;
	string inputFilePath = "Ass1Input.txt";
	inputFile.open(inputFilePath);
	string currentText;
	vector<vector<string>> puzzlePieces;
	string currentNumber;

	const int RowSize = 8;
	const int ColumnSize = 8;
	const char delimiter = ' ';

	for (int i = 0; i < RowSize; i++)
	{
		vector<string> row;

		getline(inputFile, currentText);

		stringstream currentLine(currentText);

		while (getline(currentLine, currentNumber, delimiter))
			row.push_back(currentNumber);

		puzzlePieces.push_back(row);
	}

	return puzzlePieces;
}

// fitness evaluation is adapted from professor to return a fitness score
int CountRowMismatches(vector<string> first, vector<string> second);
int CountColumnMismatches(vector<string> first, vector<string> second);

int CountRowMismatches(vector<string> first, vector<string> second)
{
	int numberOfMismatches = 0;

	for (int i = 0; i < first.size(); i++)
		if (first[i].at(2) != second[i].at(0))
			numberOfMismatches++;

	return numberOfMismatches;
}

int CountColumnMismatches(vector<string> first, vector<string> second)
{
	int numberOfMismatches = 0;

	for (int i = 0; i < first.size(); i++)
		if (first[i].at(1) != second[i].at(3))
			numberOfMismatches++;

	return numberOfMismatches;
}

int FitnessEvaluation(const vector<vector<string>>& puzzle) {
	const int RowSize = 8;
	const int ColumnSize = 8;

	int numberOfMismatches = 0;

	for (int i = 0; i < RowSize - 1; i++) {
		numberOfMismatches += CountRowMismatches(puzzle[i], puzzle[i + 1]);

	}

	for (int i = 0; i < ColumnSize - 1; i++)
	{
		vector<string> firstColumn;
		vector<string> secondColumn;

		for (int j = 0; j < RowSize; j++)
		{
			firstColumn.push_back(puzzle[j][i]);
			secondColumn.push_back(puzzle[j][i + 1]);
		}

		numberOfMismatches += CountColumnMismatches(firstColumn, secondColumn);
	}

	//inputFile.close();

	return numberOfMismatches;
}


// A Tournament Selection method to choose the best individual out of a certian tournament size picked from population
void tournamentSelection(vector<int>& parents, int populationSize, int tournamentSize, int parentSize, const vector<Individual>& population) {
	// Perform the tournament until we have selected enough parents
	for (int i = 0;i < parentSize; i++) {
		int bestIndex = -1;
		int bestFitness = INT_MAX;  // Since we want to minimize mismatches, lower fitness is better

		// Select individuals for the tournament randomly
		for (int i = 0; i < tournamentSize; i++) {
			int candidateIndex = std::rand() % populationSize;
			int candidateFitness = FitnessEvaluation(population[candidateIndex].puzzle);
			if (candidateFitness < bestFitness) {
				bestFitness = candidateFitness;
				bestIndex = candidateIndex;
			}
		}

		// Add the best individual from the tournament to the selected parents
		parents[i] = bestIndex;
	}
}

// Edge Recombination crossover Function

// Function to get the neighbors of a given piece in a puzzle which used in the edge recombination crossover function
vector<pair<int, string>> getNeighbors(const vector<vector<string>>& puzzle, int pieceIndex) {
	vector<pair<int, string>> neighbors;
	int row = pieceIndex / 8;
	int col = pieceIndex % 8;

	if (row > 0)
	{
		pair<int, string> tempPair = make_pair((row - 1) * 8 + col, puzzle[row - 1][col]);
		neighbors.push_back(tempPair);  // Above
	}

	if (row < 7) {
		pair<int, string> tempPair = make_pair((row + 1) * 8 + col, puzzle[row + 1][col]);
		neighbors.push_back(tempPair);;  // Below
	}
	if (col > 0) {
		pair<int, string> tempPair = make_pair(row * 8 + (col - 1), puzzle[row][col - 1]);
		neighbors.push_back(tempPair);  // Left
	}
	if (col < 7) {
		pair<int, string> tempPair = make_pair(row * 8 + (col + 1), puzzle[row][col + 1]);
		neighbors.push_back(tempPair);  // Right
	}

	return neighbors;
}


// Function to find common elements between two vectors of pairs which used in edge recombination crossover function
int getNumCommonElements(const std::vector<std::pair<int, std::string>>& vec1, const std::vector<std::pair<int, std::string>>& vec2) {
	std::vector<std::string> commonElements;

	// Iterate over the first vector
	for (const auto& pair1 : vec1) {
		// Check if the second element of pair1 is in the second vector
		for (const auto& pair2 : vec2) {
			if (pair1.second == pair2.second) { // Compare the strings
				commonElements.push_back(pair1.second);  // Add the common element to the result vector
				break; // Break to avoid adding duplicates
			}
		}
	}

	return commonElements.size();
}

// Function to perform Edge Recombination Crossover 
vector<vector<string>> edgeRecombinationCrossover(const vector<vector<string>>& parent1, const vector<vector<string>>& parent2) {
	unordered_map<int, vector<pair<int, std::string>>> adjacencyTable;// this is the table that used the location as the key to access the neighbors of a position of both parents


	//built adjacency Table & get tile list
	for (int i = 0; i < 64; ++i) {
		int row = i / 8;
		int col = i % 8;
		vector<pair<int, string>> tempNeighbors1 = getNeighbors(parent1, i);// store neighbors in index rather than the string content
		vector<pair<int, string>> tempNeighbors2 = getNeighbors(parent2, i);
		adjacencyTable[i].insert(adjacencyTable[i].end(), tempNeighbors1.begin(), tempNeighbors1.end());
		adjacencyTable[i].insert(adjacencyTable[i].end(), tempNeighbors2.begin(), tempNeighbors2.end());
	}

	vector<vector<string>> offspring(8, vector<string>(8, "-1"));//the returned offspring
	unordered_set<int> chosenPieces; // track already chosen pieces by position
	// Start with a random piece
	int currentPieceIndex = std::rand() % 64; // 0 to 63
	int row = currentPieceIndex / 8;
	int col = currentPieceIndex % 8;
	int pickFlag = std::rand() % 2;// 0 to 1 used for picking the first position to start with which either from parent 1 or parent 2
	pair<int, string> currentPiece;// track the current piece
	//choose the first tile from random positon from either parent1 or parent2
	if (pickFlag == 0) {
		offspring[row][col] = parent1[row][col];
	}
	else {
		offspring[row][col] = parent2[row][col];
	}

	chosenPieces.insert(currentPieceIndex);

	// Track empty positions
	vector<pair<int, int>> emptyPositions;
	for (int i = 0; i < 8; ++i) {
		for (int j = 0; j < 8; ++j) {
			if (offspring[i][j] == "-1") {
				emptyPositions.push_back({ i, j });
			}
		}
	}

	//main loop/
	for (int placedPieces = 1; placedPieces < 64; ++placedPieces) {
		vector<pair<int, string>> currentNeighbors;// neighbors of currentPiece
		// Get the neighbors based on the first chosen piece
		row = currentPieceIndex / 8;
		col = currentPieceIndex % 8;
		currentNeighbors = adjacencyTable[currentPieceIndex];

		// Remove neighbors that have already been chosen based on their position in chosenPieces
		currentNeighbors.erase(
			std::remove_if(
				currentNeighbors.begin(),
				currentNeighbors.end(),
				[&](const pair<int, string>& neighbor) {
					return chosenPieces.find(neighbor.first) != chosenPieces.end();
				}
			),
			currentNeighbors.end()
		);


		// select the neighbor with common edge, if no common edge, choose the one with the most number of edges, choose randomly if tie
		// if there is a neighbor or some neighbors
		if (!currentNeighbors.empty()) {
			int minSize = 100;
			vector<pair<int, string>> candidates1; // for most common edge
			vector < pair<int, string>> candidates2; // for most size
			for (int i = 0; i < currentNeighbors.size(); i++) {
				vector<pair<int, string>> neighbors = adjacencyTable[currentNeighbors[i].first];// neighbors of chosen currentNeighbor

				// Remove neighbors that have already been chosen based on their position in chosenPieces
				currentNeighbors.erase(
					std::remove_if(
						currentNeighbors.begin(),
						currentNeighbors.end(),
						[&](const pair<int, string>& neighbor) {
							return chosenPieces.find(neighbor.first) != chosenPieces.end();
						}
					),
					currentNeighbors.end()
				);
				
				for (int i = 0; i < currentNeighbors.size();i++) {
					// find candidates for common edge
					int numCommonElement = getNumCommonElements(currentNeighbors, neighbors);
					if (numCommonElement > 0) {
						candidates1.push_back(currentNeighbors[i]);
					}

					//find candidates for min size
					int tempSize = neighbors.size();
					if (tempSize < minSize) {
						minSize = tempSize;
						candidates2.clear();
						candidates2.push_back(currentNeighbors[i]);
					}
					if (tempSize == minSize) { //same to max
						candidates2.push_back(currentNeighbors[i]);
					}
				}
			}
			//if no common edge
			if (!candidates1.empty()) {
				// Randomly select one from the candidates1
				//std::srand(static_cast<unsigned int>(std::time(0)));  // Seed for randomness
				int randomIndex = std::rand() % candidates1.size();
				currentPiece = candidates1[randomIndex];  // select the next piece from neighbors with common edges
			}
			// then choose the one with the most edge
			else {
				// Randomly select one from the candidates2
				//std::srand(static_cast<unsigned int>(std::time(0)));  // Seed for randomness
				int randomIndex = std::rand() % candidates2.size();
				currentPiece = candidates2[randomIndex];  // select the next piece from neighbors with common edges
			}

			// Update the offspring with the chosen piece
			currentPieceIndex = currentPiece.first;
			row = currentPieceIndex / 8;
			col = currentPieceIndex % 8;
			offspring[row][col] = currentPiece.second;
			chosenPieces.insert(currentPiece.first); // Mark this piece as chosen
			//remove empty position for the one that is chosed
			emptyPositions.erase(std::remove_if(emptyPositions.begin(), emptyPositions.end(),
				[&](const pair<int, int>& pos) {
					return pos.first == row && pos.second == col;
				}), emptyPositions.end());
		}
		// if there is no neighbors for this positon, jump to a random positon, the tile is choosed randomly from either parent1 or parent2 to ensure the offspring is filled fully
		else {
			if (!emptyPositions.empty()) {
				// Jump to a random empty position
				int randomEmptyIndex = std::rand() % emptyPositions.size();
				row = emptyPositions[randomEmptyIndex].first;
				col = emptyPositions[randomEmptyIndex].second;

				// First, try to take the piece from the same position in the parents
				string selectedPiece;
				bool pieceSelected = false;

				// Generate a random integer, either 0 or 1, use to change the order of checking two parents
				int randomInt = std::rand() % 2;
				if (randomInt == 0) {
					selectedPiece = parent1[row][col];
				}
				else {
					selectedPiece = parent2[row][col];
				}

				// Place the selected piece in the offspring and mark it as chosen
				offspring[row][col] = selectedPiece;
				chosenPieces.insert(row * 8 + col);

				// Remove this position from the emptyPositions list
				emptyPositions.erase(emptyPositions.begin() + randomEmptyIndex);
			}
		}
	}
	return offspring;
}

// for testing usage
int countCommonElements(const vector<vector<string>>& matrix1, const vector<vector<string>>& matrix2) {
	unordered_set<string> uniqueElements;

	// Insert all elements of the first matrix into the set
	for (const auto& row : matrix1) {
		for (const auto& element : row) {
			uniqueElements.insert(element);
		}
	}

	int commonCount = 0;

	// Check elements of the second matrix against the set
	for (const auto& row : matrix2) {
		for (const auto& element : row) {
			if (uniqueElements.find(element) != uniqueElements.end()) {
				commonCount++;
				uniqueElements.erase(element);  // Remove to count each element only once
			}
		}
	}

	return commonCount;
}

// Function to perform swap mutation on a given puzzle matrix, swaping two random distinct position
void swapMutation(vector<vector<string>>& puzzle) {
	// Create a random number generator
	std::random_device rd;   // Non-deterministic random number generator
	std::mt19937 gen(rd());  // Mersenne Twister engine seeded with rd()

	// Uniform distribution for selecting rows and columns within the 8x8 puzzle matrix
	std::uniform_int_distribution<> dist(0, 7);

	// Randomly select two distinct positions in the puzzle matrix
	int row1, col1, row2, col2;

	do {
		row1 = dist(gen);
		col1 = dist(gen);
		row2 = dist(gen);
		col2 = dist(gen);
	} while (row1 == row2 && col1 == col2);  // Ensure that the two positions are distinct

	// Swap the elements at the two positions
	std::swap(puzzle[row1][col1], puzzle[row2][col2]);
}

// Function to rotate a tile 90 degrees clockwise, used in the roteteTwoTiles mutation function
string rotateTile90Degrees(const string& tile) {
	// Assuming the tile is a string of 4 digits representing sides [N, E, S, W]
	// A 90-degree rotation will make it [W, N, E, S]
	if (tile.size() != 4) {
		cerr << "Invalid tile format: " << tile << endl;
		return tile;  // Return the tile unchanged if it's not 4 digits
	}

	// Rotate 90 degrees clockwise: abcd �� dabc
	return string{ tile[3], tile[0], tile[1], tile[2] };
}

// Function to perform a mutation that rotates two randomly chosen tiles  90 degrees clockwise
void rotateTwoTiles(vector<vector<string>>& puzzle) {
	// Ensure we are working with an 8x8 puzzle
	if (puzzle.size() != 8 || puzzle[0].size() != 8) {
		cerr << "Puzzle is not 8x8!" << endl;
		return;
	}

	// Randomly select two distinct positions
	int tile1Row = std::rand() % 8;
	int tile1Col = std::rand() % 8;
	int tile2Row, tile2Col;

	// Ensure the second tile is different from the first
	do {
		tile2Row = std::rand() % 8;
		tile2Col = std::rand() % 8;
	} while (tile1Row == tile2Row && tile1Col == tile2Col);


	// Rotate both tiles by 90 degrees
	puzzle[tile1Row][tile1Col] = rotateTile90Degrees(puzzle[tile1Row][tile1Col]);
	puzzle[tile2Row][tile2Col] = rotateTile90Degrees(puzzle[tile2Row][tile2Col]);
}

// Function to generate a random puzzle as a grid (2D vector) of strings
vector<vector<string>> generateRandomPuzzle() {
	int gridSize = 8;//generate a grid with 8*8
	// Initialize a 2D vector of size gridSize x gridSize
	vector<vector<string>> puzzle(gridSize, vector<string>(gridSize));

	// Fill the grid with random strings, each tile will have 4 random strings (0 to 6)
	for (int row = 0; row < gridSize; ++row) {
		for (int col = 0; col < gridSize; ++col) {
			string tile = "";

			// Generate a string of 4 random characters, each between '0' and '6'
			for (int i = 0; i < 4; ++i) {
				char randomChar = '0' + rand() % 7;  // Random character between '0' and '6'
				tile += randomChar;
			}

			puzzle[row][col] = tile;  // Assign the generated tile to the puzzle
		}
	}

	return puzzle;
}

// survivor selection
// Hybrid selection based on both age and fitness, in the first half of the generations, the age is not affecting, in the second half of the generations, older individual has some advantage
void hybridSurvivorSelection(vector<Individual>& population, int currentGeneration, int maxGenerations, int populationSize) {
	// Parameters for hybrid selection
	double transition_point = 0.5;  // Start transitioning from low age to high age after 50% of generations
	double age_weight = 0.0;         // No weight on age in the early generations
	double fitness_weight = 1.0;     // Full weight on fitness in the early generations

	// Calculate current generation's age vs fitness weighting
	double generation_ratio = static_cast<double>(currentGeneration) / maxGenerations;

	//elitism selection logic
	// Sort population based on a combined score of age and fitness
	std::sort(population.begin(), population.end(), [&](const Individual& a, const Individual& b) {
		double a_score = age_weight * (1.0 / a.age) + fitness_weight * a.fitness;  // Lower age is better
		double b_score = age_weight * (1.0 / b.age) + fitness_weight * b.fitness;  // Lower age is better
		return a_score < b_score;  // Sorting in ascending order
		});

	// Keep only the best individuals based on the desired population size
	if (population.size() > populationSize) {
		population.resize(populationSize);
	}
}

// Update age of each individual in the population
void updateAge(vector<Individual>& population) {
	for (auto& individual : population) {
		individual.age += 1;  // Increment age of each individual
	}
}

// function to generate the best individual amoung the current population
Individual generateBestIndividual(const vector<Individual> population) {
	int bestFitness = numeric_limits<int>::max();  // Start with the lowest possible value
	Individual bestIndividual;

	// Iterate through each individual in the population
	for (const auto& individual : population) {
		// Calculate the fitness of the current individual
		double fitness = individual.fitness;

		// If the current fitness is better than the best fitness found so far, update
		if (fitness < bestFitness) {
			bestFitness = fitness;
			bestIndividual = individual;  // Store the current individual as the best one
		}
	}

	return bestIndividual;  // Return the individual with the best fitness
}

// function to output a individual to a txt file
void outputIndividualToFile(const Individual& ind, const string& filename) {
	// Open a file for writing
	ofstream outFile(filename);

	// Check if the file is open
	if (!outFile) {
		cerr << "Error opening file: " << filename << endl;
		return;
	}
	outFile << "Qiaoyu Xu 40292010";
	outFile << endl;

	// Write puzzle to file
	for (size_t i = 0; i < ind.puzzle.size(); ++i) {
		for (size_t j = 0; j < ind.puzzle[i].size(); ++j) {
			outFile << ind.puzzle[i][j];  // Write each element in the row
			if (j != ind.puzzle[i].size() - 1) {
				outFile << " ";  // Add space only between elements, not after the last one
			}
		}
		if (i != ind.puzzle.size() - 1) {
			outFile << endl;  // New line for the next row except the last one
		}
	}

	// Close the file
	outFile.close();
}

int main()
{
	// Seed the random number generator
	std::srand(static_cast<unsigned int>(std::time(0)));
	
	// Ask the user how many files to generate
	int populationSize;
	std::cout << "Enter the number of populaton size: ";
	std::cin >> populationSize;

	// Ask the user how many generation
	int numberOfGeneration;
	std::cout << "Enter the number of generation to you want: ";
	std::cin >> numberOfGeneration;

	int parentSize = populationSize*0.5;
	int tournamentSize = 0;//for tournament selection
	double tempTournamentSize = populationSize * 0.1;// tournament size is 10% of population size initially
	tournamentSize = static_cast<int>(std::ceil(tempTournamentSize));
	double mutationProbability = 0.3;

	// Initialize population
	vector<Individual> population(populationSize);
	// get the ass1input into the population
	population[0].puzzle = readFile();
	population[0].fitness = FitnessEvaluation(population[0].puzzle);
	population[0].age = 0;  // Everyone starts with age 0

	// Generate rest of population after reading the Ass1Input.txt
	for (int i = 1; i < populationSize; ++i) {
		population[i].puzzle = generateRandomPuzzle();
		population[i].fitness = FitnessEvaluation(population[i].puzzle);
		population[i].age = 0;  // Everyone starts with age 0
	}

	//initialize the bestIndividual
	Individual bestInvididual;
	bestInvididual.fitness = 99999;
	bestInvididual.age = 0;

	// Evolve population over generations
	for (int generation = 0; generation < numberOfGeneration; ++generation) {
		std::cout << "Generation " << generation + 1 << std::endl;

		// Parent selection (Tournament or other selection method)
		vector<int> parents(parentSize);
		tournamentSelection(parents, populationSize, tournamentSize, parentSize, population);
		double generation_ratio = static_cast<double>(generation) / numberOfGeneration;
		mutationProbability = 0.3 - 0.2 * generation_ratio;

		// Crossover and mutation to generate offspring
		for (int i = 0; i < parentSize / 2; ++i) {
			vector<vector<string>> parent1 = population[parents[2 * i]].puzzle;
			vector<vector<string>> parent2 = population[parents[2 * i + 1]].puzzle;

			// Perform crossover
			vector<vector<string>> offspring = edgeRecombinationCrossover(parent1, parent2);

			// Apply mutation with a certain probability (e.g., 0.1 or 10%)
			double randProb = static_cast<double>(rand()) / RAND_MAX;
			if (randProb < mutationProbability) {
				// Randomly choose which mutation to apply
				if (rand() % 2 == 0) {
					swapMutation(offspring);  // Apply swap mutation
				}
				else {
					rotateTwoTiles(offspring); // Apply rotate mutation
				}
			}
			// Add the offspring to the population
			population.push_back({ offspring, FitnessEvaluation(offspring), 0 });  // Age is 0 for new offspring
		}

		// Update the age of all individuals in the population
		updateAge(population);

		// Apply hybrid survivor selection
		hybridSurvivorSelection(population, generation, numberOfGeneration, populationSize);
		Individual candidateBestInvididual = generateBestIndividual(population);
		if (candidateBestInvididual.fitness < bestInvididual.fitness) {
			bestInvididual = candidateBestInvididual;
		}
	}

	cout << "the fitness of the best individual is " << bestInvididual.fitness;
	//output the best individual to txt file
	outputIndividualToFile(bestInvididual, "Ass1Output.txt");


	return 0;
}