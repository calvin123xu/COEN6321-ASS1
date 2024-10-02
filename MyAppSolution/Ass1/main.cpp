#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <iomanip>    // For formatting output
#include <cstdlib>    // For rand() and srand()
#include <ctime>      // For time()
#include <unordered_map>
#include <algorithm>  // For std::min_element and std::unique
#include <utility>  // For std::pair
#include <unordered_set>
#include <functional>

using namespace std;

// Custom hash function for std::pair<int, std::string>
struct pair_hash {
	template <class T1, class T2>
	std::size_t operator () (const std::pair<T1, T2>& pair) const {
		auto hash1 = std::hash<T1>{}(pair.first);     // Hash the first element
		auto hash2 = std::hash<T2>{}(pair.second);    // Hash the second element
		return hash1 ^ (hash2 << 1);                  // Combine the two hashes
	}
};

// Custom equality comparator for std::pair<int, std::string>
struct pair_equal {
	bool operator() (const std::pair<int, std::string>& lhs, const std::pair<int, std::string>& rhs) const {
		return lhs.first == rhs.first && lhs.second == rhs.second; // Compare both elements
	}
};

// Function to generate random 4-digit numbers between 0 and 6 in to a txt file
void generateFile(int fileNumber) {
	// Create the file name "outputy.txt"
	std::string directory = "outputs/";
	std::string fileName = directory + "output" + std::to_string(fileNumber) + ".txt";

	// Open a file to write the output
	std::ofstream outFile(fileName);
	if (!outFile) {
		std::cerr << "Error opening file: " << fileName << std::endl;
		return;
	}

	// Create a 2D array to store the random numbers
	int numbers[8][8];

	// First, generate the sequence of 8 rows and 8 columns of 4-digit numbers
	for (int row = 0; row < 8; ++row) {
		for (int col = 0; col < 8; ++col) {
			int num = 0;
			for (int i = 0; i < 4; ++i) {
				num = num * 10 + (std::rand() % 7);  // % 7 generates a number from 0 to 6
			}
			numbers[row][col] = num;  // Store the number in the array
		}
	}

	// Now, randomly pick some numbers and modify them
	int modifications = std::rand() % 10 + 5;  // Randomly choose between 5 and 14 numbers to modify
	for (int i = 0; i < modifications; ++i) {
		int randRow = std::rand() % 8;
		int randCol = std::rand() % 8;
		int newNum = 0;
		for (int j = 0; j < 4; ++j) {
			newNum = newNum * 10 + (std::rand() % 7);  // Generate a new random 4-digit number
		}
		numbers[randRow][randCol] = newNum;  // Modify the selected number
	}

	// Write the modified sequence to the file
	for (int row = 0; row < 8; ++row) {
		for (int col = 0; col < 7; ++col) {
			outFile << std::setw(4) << std::setfill('0') << numbers[row][col] << " ";
		}
		outFile << std::setw(4) << std::setfill('0') << numbers[row][7];
		outFile << std::endl;  // Newline at the end of each row
	}

	// Close the file
	outFile.close();
}


// fitness evaluation
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

int FitnessEvaluation(int fileNumber) {
	//fitness evaluation 
	ifstream inputFile;
	std::string directory = "outputs/";
	string inputFilePath = directory + "output" + std::to_string(fileNumber) + ".txt";
	inputFile.open(inputFilePath);
	string currentText;
	//getline(inputFile, currentText);////////////////////////////modification

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

	int numberOfMismatches = 0;

	for (int i = 0; i < RowSize - 1; i++) {
		numberOfMismatches += CountRowMismatches(puzzlePieces[i], puzzlePieces[i + 1]);

	}

	for (int i = 0; i < ColumnSize - 1; i++)
	{
		vector<string> firstColumn;
		vector<string> secondColumn;

		for (int j = 0; j < RowSize; j++)
		{
			firstColumn.push_back(puzzlePieces[j][i]);
			secondColumn.push_back(puzzlePieces[j][i + 1]);
		}

		numberOfMismatches += CountColumnMismatches(firstColumn, secondColumn);
	}

	inputFile.close();

	return numberOfMismatches;
}

// Tournament Selection Function giving tournamentSize and parentSize, later to implement dynamically change tournamentSize and parentSize according to generation and diversity
void tournamentSelection(vector<int>& parents, int populationSize, int tournamentSize, int parentSize) {
	// Perform the tournament until we have selected enough parents
	for (int i = 0;i < parentSize; i++) {
		int bestIndex = -1;
		int bestFitness = INT_MAX;  // Since we want to minimize mismatches, lower fitness is better

		// Select individuals for the tournament randomly
		for (int i = 0; i < tournamentSize; i++) {
			int candidateIndex = std::rand() % populationSize;
			int candidateFitness = FitnessEvaluation(candidateIndex);
			if (candidateFitness < bestFitness) {
				bestFitness = candidateFitness;
				bestIndex = candidateIndex;
			}
		}

		// Add the best individual from the tournament to the selected parents
		parents[i] = bestIndex;
	}
}

// Helper function to read the puzzle file into a 2D vector of strings
vector<vector<string>> readPuzzle(int fileNumber) {
	ifstream inputFile;
	std::string directory = "outputs/";
	string inputFilePath = directory + "output" + std::to_string(fileNumber) + ".txt";
	inputFile.open(inputFilePath);
	string currentText;

	vector<vector<string>> puzzlePieces;
	string currentNumber;
	const char delimiter = ' ';

	for (int i = 0; i < 8; i++) {
		vector<string> row;
		getline(inputFile, currentText);
		stringstream currentLine(currentText);
		while (getline(currentLine, currentNumber, delimiter))
			row.push_back(currentNumber);
		puzzlePieces.push_back(row);
	}

	inputFile.close();
	return puzzlePieces;
}

// ERX Functions

// Function to get the neighbors of a given piece in a puzzle
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



// for test purpose
void printVector(const std::vector<std::vector<std::string>>& vec) {
	for (const auto& row : vec) {         // Iterate over each row
		for (const auto& element : row) {  // Iterate over each element in the row
			std::cout << element << " ";    // Print the element followed by a space
		}
		std::cout << std::endl;             // Print a newline after each row
	}
}

// Function to find common elements between two vectors of pairs
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
	unordered_map<std::pair<int, std::string>,vector<pair<int, std::string>>, pair_hash, pair_equal> adjacencyTable;
	vector<string> tileList;// list all the possible tiles in both parent
	//built adjacency Table & get tile list
	for (int i = 0; i < 64; ++i) {
		int row = i / 8;
		int col = i % 8;
		vector<pair<int, string>> tempNeighbors1 = getNeighbors(parent1, i);// store neighbors in index rather than the string content
		vector<pair<int, string>> tempNeighbors2 = getNeighbors(parent2, i);
		pair<int, string > tempPair1 = make_pair(i, parent1[row][col]);
		pair<int, string > tempPair2 = make_pair(i, parent2[row][col]);
		adjacencyTable[tempPair1].insert(adjacencyTable[tempPair1].end(), tempNeighbors1.begin(), tempNeighbors1.end());
		adjacencyTable[tempPair2].insert(adjacencyTable[tempPair2].end(), tempNeighbors2.begin(), tempNeighbors2.end());
		tileList.push_back(parent1[row][col]);
		tileList.push_back(parent2[row][col]);
	}
	
	// Eliminate duplicate elements in tileList
	std::sort(tileList.begin(), tileList.end());
	auto it = std::unique(tileList.begin(), tileList.end());
	tileList.erase(it, tileList.end());
	

	vector<vector<string>> offspring(8, vector<string>(8, "-1"));//the returned offspring
	unordered_set<string> chosenPieces; // track already chosen pieces
	// Start with a random piece
	int currentPieceIndex = std::rand() % 64; // 0 to 63
	int row = currentPieceIndex / 8;
	int col = currentPieceIndex % 8;
	int pickFlag = std::rand() % 2;// 0 to 1
	pair<int,string> currentPiece;// track the current piece

	if (pickFlag == 0) {
		offspring[row][col] = parent1[row][col];
	}
	else {
		offspring[row][col] = parent2[row][col];
	}

	chosenPieces.insert(offspring[row][col]);

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
	int i = 1;//tracking number of time in main loop
	//for (int placedPieces = 1; placedPieces < 64; ++placedPieces) {
	while (chosenPieces.size() < 64){
		vector<pair<int, string>> currentNeighbors;// neighbors of currentPiece

		if (i == 1) {// for first iteration
			if (pickFlag == 0) {
				currentNeighbors = adjacencyTable[make_pair(currentPieceIndex, parent1[row][col])];
			}
			else {
				currentNeighbors = adjacencyTable[make_pair(currentPieceIndex, parent1[row][col])];
			}
		}
		else {
			// Get the neighbors based on the first chosen piece
			row = currentPieceIndex / 8;
			col = currentPieceIndex % 8;
			currentNeighbors = adjacencyTable[make_pair(currentPieceIndex, offspring[row][col])];
		}

		i++;

		// Remove neighbors that have already been chosen in the offspring
		currentNeighbors.erase(std::remove_if(currentNeighbors.begin(), currentNeighbors.end(),
			[&](const pair<int, string>& neighbor) {
				return chosenPieces.find(neighbor.second) != chosenPieces.end();
			}), currentNeighbors.end());

		// select the neighbor with the most edge, choose randomly if tie
		if (!currentNeighbors.empty()) {
			int maxCommonElement = 0;
			int maxSize = -1;
			vector<pair<int,string>> candidates1; // for most common edge
			vector < pair<int,string>> candidates2; // for most size
			for (int i = 0; i < currentNeighbors.size(); i++) {
				vector<pair<int, string>> neighbors = adjacencyTable[currentNeighbors[i]];// neighbors of chosen currentNeighbor

				//erase neighbors what has already bean chosed in offspring
				// Remove neighbors that have already been chosen in the offspring
				neighbors.erase(std::remove_if(neighbors.begin(), neighbors.end(),
					[&](const pair<int, string>& neighbor) {
						return chosenPieces.find(neighbor.second) != chosenPieces.end();
					}), neighbors.end());

				for (int i = 0; i < currentNeighbors.size();i++) {
					// find candidates for most common edges
					int numCommonElement = getNumCommonElements(currentNeighbors, neighbors);
					if (numCommonElement > maxCommonElement) {
						maxCommonElement = numCommonElement;
						candidates1.clear();
						candidates1.push_back(currentNeighbors[i]);
					}
					if (numCommonElement = maxCommonElement) { //same to max
						candidates1.push_back(currentNeighbors[i]);
					}
					//find candidates for most size
					int tempSize = neighbors.size();
					if (tempSize > maxSize) {
						maxSize = tempSize;
						candidates2.clear();
						candidates2.push_back(currentNeighbors[i]);
					}
					if (numCommonElement = maxCommonElement) { //same to max
						candidates2.push_back(currentNeighbors[i]);
					}
				}
			}
			if (!candidates1.empty()) {
				// Randomly select one from the candidates1
				std::srand(static_cast<unsigned int>(std::time(0)));  // Seed for randomness
				int randomIndex = std::rand() % candidates1.size();
				currentPiece = candidates1[randomIndex];  // select the next piece from neighbors with common edges
			}
			else {
				// Randomly select one from the candidates2
				std::srand(static_cast<unsigned int>(std::time(0)));  // Seed for randomness
				int randomIndex = std::rand() % candidates2.size();
				currentPiece = candidates2[randomIndex];  // select the next piece from neighbors with common edges
			}
			// Update the offspring with the chosen piece
			currentPieceIndex = currentPiece.first;
			row = currentPieceIndex / 8;
			col = currentPieceIndex % 8;
			offspring[row][col] = currentPiece.second;
			chosenPieces.insert(currentPiece.second); // Mark this piece as chosen
			emptyPositions.erase(std::remove_if(emptyPositions.begin(), emptyPositions.end(),
				[&](const pair<int, int>& pos) {
					return pos.first == row && pos.second == col;
				}), emptyPositions.end());
		}
		else {
			if (!emptyPositions.empty()) {
				// Jump to a random empty position
				int randomEmptyIndex = std::rand() % emptyPositions.size();
				row = emptyPositions[randomEmptyIndex].first;
				col = emptyPositions[randomEmptyIndex].second;

				// First, try to take the piece from the same position in the parents
				string selectedPiece;
				bool pieceSelected = false;

				// Seed the random number generator with the current time
				std::srand(std::time(0));
				// Generate a random integer, either 0 or 1, use to change the order of checking two parents
				int randomInt = std::rand() % 2;
				if (randomInt == 0) {
					// Check if the piece at this position in parent1 is available
					if (chosenPieces.find(parent1[row][col]) == chosenPieces.end()) {
						selectedPiece = parent1[row][col];
						pieceSelected = true;
					}
					// Check if the piece at this position in parent2 is available
					else if (chosenPieces.find(parent2[row][col]) == chosenPieces.end()) {
						selectedPiece = parent2[row][col];
						pieceSelected = true;
					}
				}
				else {
					// Check if the piece at this position in parent2 is available
					if (chosenPieces.find(parent2[row][col]) == chosenPieces.end()) {
						selectedPiece = parent2[row][col];
						pieceSelected = true;
					}
					// Check if the piece at this position in parent1 is available
					else if (chosenPieces.find(parent1[row][col]) == chosenPieces.end()) {
						selectedPiece = parent1[row][col];
						pieceSelected = true;
					}
				}
				

				// If neither piece from the parents at the same position is available, choose randomly from tileList
				if (!pieceSelected) {
					do {
						selectedPiece = tileList[std::rand() % tileList.size()];
					} while (chosenPieces.find(selectedPiece) != chosenPieces.end());
				}

				// Place the selected piece in the offspring and mark it as chosen
				offspring[row][col] = selectedPiece;
				chosenPieces.insert(selectedPiece);

				// Remove this position from the emptyPositions list
				emptyPositions.erase(emptyPositions.begin() + randomEmptyIndex);
			}
		}
		
		
		printVector(offspring);
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

	// Generate x files with names output0.txt, output1.txt, ..., outputx-1.txt
	for (int i = 0; i < populationSize; ++i) {
		generateFile(i);  // Generate and modify the i-th file
	}

	//tournament selection to pick a pool of parents, the tournament size of 5% of the population which round up to a interger, the parentSize is currently set to 50
	int tournamentSize = 0;
	double tempTournamentSize = populationSize * 0.1;
	tournamentSize = static_cast<int>(std::ceil(tempTournamentSize));
	vector<int> parents(populationSize);
	tournamentSelection(parents, populationSize, tournamentSize, 50);

	//offspring generation
	// edge recombination crossover
	//pick two parents in the pool by ranking selection
	 // Read the two selected parent puzzles
	vector<vector<string>> parent1 = readPuzzle(0);
	vector<vector<string>> parent2 = readPuzzle(1);

	// Perform Edge Recombination Crossover to generate offspring
	vector<vector<string>> offspring = edgeRecombinationCrossover(parent1, parent2);

	// Output the offspring puzzle
	cout << "Offspring Puzzle:" << endl;
	for (const auto& row : offspring) {
		for (const auto& piece : row) {
			cout << piece << " ";
		}
		cout << endl;
	}
	cout << "Number of elements in common: " << to_string(countCommonElements(parent1, parent2)) << endl;

	return 0;
}