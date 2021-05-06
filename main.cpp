#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

using namespace std;

#define DEBUGVAR(x) do { std::cerr << #x << ": " << x << std::endl; } while (0)

/**
 * Auto-generated code below aims at helping you parse
 * the standard input according to the problem statement.
 **/

int main()
{
    int numberOfCells; // 37
    cin >> numberOfCells; cin.ignore();
    for (int i = 0; i < numberOfCells; i++) {
        int index; // 0 is the center cell, the next cells spiral outwards
        int richness; // 0 if the cell is unusable, 1-3 for usable cells
        int neigh0; // the index of the neighbouring cell for each direction
        int neigh1;
        int neigh2;
        int neigh3;
        int neigh4;
        int neigh5;
        cin >> index >> richness >> neigh0 >> neigh1 >> neigh2 >> neigh3 >> neigh4 >> neigh5; cin.ignore();
    }

    // game loop
    while (1) {
        int day; // the game lasts 24 days: 0-23
        cin >> day; cin.ignore();
        int nutrients; // the base score you gain from the next COMPLETE action
        cin >> nutrients; cin.ignore();
        int sun; // your sun points
        int score; // your current score
        cin >> sun >> score; cin.ignore();
        int oppSun; // opponent's sun points
        int oppScore; // opponent's score
        bool oppIsWaiting; // whether your opponent is asleep until the next day
        cin >> oppSun >> oppScore >> oppIsWaiting; cin.ignore();
        int numberOfTrees; // the current amount of trees
        cin >> numberOfTrees; cin.ignore();
        pair<int, int> chosenTree(0, 0);
        for (int i = 0; i < numberOfTrees; i++) {
            int cellIndex; // location of this tree
            int size; // size of this tree: 0-3
            bool isMine; // 1 if this is your tree
            bool isDormant; // 1 if this tree is dormant
            cin >> cellIndex >> size >> isMine >> isDormant; cin.ignore();
            if (isMine && size > chosenTree.second) {
                chosenTree.first = cellIndex;
                chosenTree.second = size;
            }
        }
        int numberOfPossibleMoves;
        cin >> numberOfPossibleMoves; cin.ignore();
        DEBUGVAR(numberOfPossibleMoves);
        for (int i = 0; i < numberOfPossibleMoves; i++) {
            string possibleMove;
            getline(cin, possibleMove);
        }
        cout << "COMPLETE " << chosenTree.first << endl;
        // Write an action using cout. DON'T FORGET THE "<< endl"
        // To debug: cerr << "Debug messages..." << endl;


        // GROW cellIdx | SEED sourceIdx targetIdx | COMPLETE cellIdx | WAIT <message>
        // cout << "WAIT" << endl;
    }
}