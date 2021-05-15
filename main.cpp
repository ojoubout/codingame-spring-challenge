#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <array>
#include <set>
#include <chrono>
#include <cmath>

using namespace std;

#define DEBUGVAR(x) do { std::cerr << #x << ": " << x << std::endl; } while (0)
#define MAX_TREE 8
#define NOW chrono::high_resolution_clock::now()
#define CMLT_EXP ceil((23 - day) / 2.0)

auto cmp = [](pair<int, int> a, pair<int, int> b) { return a.second > b.second; };
auto cmp2 = [](pair<int, pair<int, int>> a, pair<int, pair<int, int>> b) { return a.second.second > b.second.second; };

int day; // the game lasts 24 days: 0-23
int sun; // your sun points

class Point {
public:
    int x;
    int y;
    int z;

    Point() {};
    Point(int x, int y, int z) : x(x), y(y), z(z) {};
};

class Tree;
// class Cell;

// Cell *cells;
// array<Tree, 37> trees;

class Cell {
public:
    static int numberOfCells;
    static array<Cell, 37> cells;

    int index;
    Point coord;
    int richness;
    int neigh[6];
    Tree *tree;

    Cell *getNeigh(int index) {
        Cell *ret;

        ret = neigh[index % 6] != -1 ? &cells[neigh[index % 6]] : NULL;

        return ret;
    }

    void print() {
        cerr << "Cell\n\t";
        DEBUGVAR(index);
        cerr << '\t';
        DEBUGVAR(richness);
        cerr << "\tneigh: " << neigh[0] << " " << neigh[1] << " " << neigh[2] << " " << neigh[3] << " "
            << neigh[4] << " " << neigh[5] << '\n';
    }
};

int Cell::numberOfCells = 0;
array<Cell, 37> Cell::cells;
int cellDistance(const Cell &c1, const Cell &c2) {
    return (abs(c1.coord.x - c2.coord.x) + abs(c1.coord.y - c2.coord.y) + abs(c1.coord.z - c2.coord.z)) / 2;
}

int canShadow(const Cell &c1, const Cell &c2) {
    if (c1.coord.x == c2.coord.x || c1.coord.y == c2.coord.y
    || c1.coord.z == c2.coord.z) {
        int dis = cellDistance(c1, c2);
        if (dis <= 3) {
            return dis;
        }
    }
    return 0;
}

class Tree {
public:
    static int numberOfTrees;
    static array<Tree, 37> trees;

    Cell *cell;
    int size;
    bool isMine;
    bool isDormant;

    Tree() : cell(NULL), size(-1) {
    }

    bool exist() const {
        return cell != NULL;
    }

    int growCost() const {
        int cost = 0;
        if (size < 3) {
            cost = size * (size + 1) + 1 + countMyTrees(size + 1);
        }
        return cost;
    }

    // pair<int, int> getNoShadowCell() {
    //     set<int> resultList;
    //     // set<int> result3;
    //     for (int i = 0; i < 6; i++) {
    //         if (cell->getNeigh(i)) {
    //             Cell *neigh1;
    //             neigh1 = cell->getNeigh(i)->getNeigh(i + 1);
    //             if (neigh1 && neigh1->tree == NULL && neigh1->richness > 0) {
    //                 resultList.insert(neigh1->index);
    //             }
    //             if (neigh1) {
    //                 Cell *neigh2, *neigh3;
    //                 neigh2 = neigh1->getNeigh(i);
    //                 neigh3 = neigh1->getNeigh(i + 1);
    //                 if (neigh2 && neigh2->tree == NULL && neigh2->richness > 0) {
    //                     resultList.insert(neigh2->index);
    //                 }
    //                 if (neigh3 && neigh3->tree == NULL && neigh3->richness > 0) {
    //                     resultList.insert(neigh3->index);
    //                 }
    //             }
    //         }
    //     }
    //     pair<int, int> result(-1, -100);
    //     for (int i : resultList) {
    //         int score = Cell::cells[i].richness * 6 - (cellDistance(*cell, Cell::cells[i]) / 3);
    //         for (const Tree &tree : trees) {
    //             // cerr << tree.cell->index << endl;
    //             if (tree.exist() && tree.isMine && canShadow(Cell::cells[i], *tree.cell)) {
    //                 // cerr << "SHADOW " << i << " " << tree.cell->index << endl;
    //                 score -= 10;
    //             }
    //         }
    //         // cerr << "Try " << i << " " << score << " " << cell->index << " " << Cell::cells[i].index << " " << cellDistance(*cell, Cell::cells[i]) << endl;
    //         if (result.first == -1 || result.second < score || (result.second == score && 
    //         cellDistance(*cell, Cell::cells[i]) < cellDistance(*cell, Cell::cells[result.first]))) {
    //             // cerr << "OK" << endl;
    //             result = make_pair(i, score);
    //         }
    //         // cerr << i << " ";
    //     }
    //     // cerr << "\n";
    //     // for (int i : result3) {
    //     //     cerr << i << " ";
    //     // }
    //     // if (result2.begin() != result2.end()) {
    //     //     return cells + (*result2.begin());
    //     // }
    //     // if (result2.begin() != result2.end()) {
    //     //     return cells + (*result2.begin());
    //     // }
    //     // return result.first != -1 ? &Cell::cells[result.first] : NULL;
    //     return result;
    // }


    static set<pair<int, pair<int, int>>, decltype(cmp2)> getSeedPlant() {
        set<pair<int, pair<int, int>>, decltype(cmp2)> result(cmp2);
        for (const Tree &tree : trees) {
            if (tree.exist() && tree.isMine) {
                pair<int, int> best(-1, -1000);
                for (int i = 0; i < 6; i++) {
                    Cell *neigh[6] = {NULL};
                    neigh[0] = tree.cell->getNeigh(i);
                    neigh[1] = neigh[0] ? neigh[0]->getNeigh(i) : NULL;
                    neigh[2] = neigh[1] ? neigh[1]->getNeigh(i) : NULL;
                    neigh[3] = neigh[0] ? neigh[0]->getNeigh(i + 1) : NULL;
                    neigh[4] = neigh[1] ? neigh[1]->getNeigh(i + 1) : NULL;
                    neigh[5] = neigh[3] ? neigh[3]->getNeigh(i + 1) : NULL;
                    for (int j = 0; j < 6; j++) {
                        if (neigh[j] && neigh[j]->tree == NULL && neigh[j]->richness > 0) {
                            int score = neigh[j]->richness - (cellDistance(*tree.cell, *neigh[j]) / 3);
                            for (const Tree &tree1 : trees) {
                                if (tree.cell != tree1.cell && tree1.exist() && tree1.isMine) {
                                    score += cellDistance(*neigh[j], *tree1.cell) * 1.5;
                                }
                                if (tree1.exist() && tree1.isMine && canShadow(*neigh[j], *tree1.cell)) {
                                    score -= 20;
                                }
                            }

                            cerr << tree.cell->index << " " << neigh[j]->index << " " << score << endl;
                            if (best.first == -1 || best.second < score) {
                                best = make_pair(neigh[j]->index, score);
                            }
                        }
                    }
                }
                result.insert(make_pair(tree.cell->index, best));
            }
        }
        return result;
    }

    static set<pair<int, int>, decltype(cmp)> getGrowSeed() {

        set<pair<int, int>, decltype(cmp)> result(cmp);
        // pair<int, int> best(-1, -100);
        for (const Tree &seed : trees) {
            if (seed.exist() && seed.isMine && seed.size < 3) {
                int score = seed.size == 0 ? 2 : seed.size - countMyTrees(seed.size + 1)*2;
                // for (const Tree &tree : trees) {
                //     if (tree.exist() && tree.isMine && canShadow(*seed.cell, *tree.cell)) {
                //         // cerr << "CAN SHADOW " << seed.cell->index << " " << tree.cell->index << endl;
                //         --score;
                //     }
                // }
                result.insert(make_pair(seed.cell->index, score));
                // if (best.first == -1 || best.second < score) {
                //     best = make_pair(seed.cell->index, score);
                // }
            }
        }
        return result;
    }

    static set<pair<int, int>, decltype(cmp)> getCompleteTree() {
        // pair<int, int> best(-1, -100);
        set<pair<int, int>, decltype(cmp)> result(cmp);

        for (const Tree &seed : trees) {
            if (seed.exist() && seed.isMine && seed.size == 3) {
                int score = -seed.cell->richness;
                for (const Tree &tree : trees) {
                    if (tree.exist() && tree.isMine && canShadow(*seed.cell, *tree.cell)) {
                        // cerr << "CAN SHADOW " << seed.cell->index << " " << tree.cell->index << endl;
                        ++score;
                    }
                }
                result.insert(make_pair(seed.cell->index, score));
            }
        }
        return result;
    }

    static int countMyTrees(int size) {
        int count = 0;
        for (const Tree &tree : trees) {
            if (tree.exist() && tree.isMine && (size == -1 || tree.size == size)) {
                ++count;
            }
        }
        return count;
    }
};


Point cellsCoordinates[37] = {
Point(0, 0, 0), Point(1, -1, 0), Point(1, 0, -1), Point(0, 1, -1), Point(-1, 1, 0),
Point(-1, 0, 1), Point(0, -1, 1), Point(2, -2, 0), Point(2, -1, -1), Point(2, 0, -2),
Point(1, 1, -2), Point(0, 2, -2), Point(-1, 2, -1), Point(-2, 2, 0), Point(-2, 1, 1),
Point(-2, 0, 2), Point(-1, -1, 2), Point(0, -2, 2), Point(1, -2, 1), Point(3, -3, 0),
Point(3, -2, -1), Point(3, -1, -2), Point(3, 0, -3), Point(2, 1, -3), Point(1, 2, -3),
Point(0, 3, -3), Point(-1, 3, -2), Point(-2, 3, -1), Point(-3, 3, 0), Point(-3, 2, 1),
Point(-3, 1, 2), Point(-3, 0, 3), Point(-2, -1, 3), Point(-1, -2, 3), Point(0, -3, 3),
Point(1, -3, 2), Point(2, -3, 1)};

int Tree::numberOfTrees = 0;
array<Tree, 37> Tree::trees;
/**
 * Auto-generated code below aims at helping you parse
 * the standard input according to the problem statement.
 **/

bool completeTrees();

namespace action {
    bool seed(Tree & tree, Cell & cell) {
        if (!tree.isDormant) {
            cout << "SEED " << tree.cell->index << " " << cell.index << endl;
            return true;
        }
        cerr << "IS DORMANT" << endl;
        return false;
    }

    bool grow(Tree & tree) {
        if (!tree.isDormant) {
            if (tree.size == 2 && ((Tree::countMyTrees(3) + 1) >= CMLT_EXP)
            && sun >= (tree.growCost() + 4)) {
                return !completeTrees();
            } else {
                cout << "GROW " << tree.cell->index << endl;
                return true;
            }
        }
        return false;
    }

    bool complete(Tree & tree) {
        if (!tree.isDormant) {
            cout << "COMPLETE " << tree.cell->index << endl;
            return true;
        }
        return false;
    }

    void wait() {
        cout << "WAIT" << endl;
    }
}


bool completeTrees() {
    cerr << "C " << Tree::countMyTrees(3) << " " << CMLT_EXP << endl;
    set<pair<int, int>, decltype(cmp)> trees = Tree::getCompleteTree();
    int i = 0;
    for (const pair<int, int> &tree : trees) {
        cerr << "CT " << tree.first << endl;
        if (action::complete(Tree::trees[tree.first])) {
            break;
        }
        ++i;
    }
    if (i == trees.size()) {
        return true;
    }
    return false;
}

double easeOutQuint(double x) {
    return 1.0 - pow(1.0 - x, 5);
}

int main()
{
    cin >> Cell::numberOfCells; cin.ignore();
    int cellsRichness[Cell::numberOfCells];
    // cells = new Cell[Cell::numberOfCells];
    for (int i = 0; i < Cell::numberOfCells; i++) {
        int index; // 0 is the center cell, the next cells spiral outwards
        cin >> index
            >> Cell::cells[index].richness
            >> Cell::cells[index].neigh[0]
            >> Cell::cells[index].neigh[1]
            >> Cell::cells[index].neigh[2]
            >> Cell::cells[index].neigh[3]
            >> Cell::cells[index].neigh[4]
            >> Cell::cells[index].neigh[5];
            cin.ignore();
        Cell::cells[index].index = index;
        Cell::cells[index].coord = cellsCoordinates[i];
        Tree::trees[index].cell = NULL;
        // Cell::cells[index].print();
        // cerr << index << "\t" << Cell::cells[index].coord.x << "\t" << Cell::cells[index].coord.y
        //     << "\t" << Cell::cells[index].coord.z << "\t|\t"
        //     << (Cell::cells[index].coord.x + Cell::cells[index].coord.y + Cell::cells[index].coord.z) << "\n";
    }

    // game loop

    while (1) {
        chrono::high_resolution_clock::time_point start = NOW;
        cin >> day;
        cin.ignore();
        int nutrients; // the base score you gain from the next COMPLETE action
        cin >> nutrients; cin.ignore();
        int score; // your current score
        cin >> sun >> score; cin.ignore();
        int oppSun; // opponent's sun points
        int oppScore; // opponent's score
        bool oppIsWaiting; // whether your opponent is asleep until the next day
        cin >> oppSun >> oppScore >> oppIsWaiting; cin.ignore();
        cin >> Tree::numberOfTrees; cin.ignore();


        for (int i = 0; i < Cell::numberOfCells; i++) {
            Tree::trees[i].cell = NULL;
            Cell::cells[i].tree = NULL;
        }

        for (int i = 0; i < Tree::numberOfTrees; i++) {
            int cellIndex; // location of this tree
            int size; // size of this tree: 0-3
            bool isMine; // 1 if this is your tree
            bool isDormant; // 1 if this tree is dormant
            cin >> cellIndex
                >> Tree::trees[cellIndex].size
                >> Tree::trees[cellIndex].isMine
                >> Tree::trees[cellIndex].isDormant;
            cin.ignore();
            Tree::trees[cellIndex].cell = &Cell::cells[cellIndex];
            Cell::cells[cellIndex].tree = &Tree::trees[cellIndex];

            // cerr << "Tree " << cellIndex << endl;
            // if (Tree::trees[cellIndex].isMine) {
            //     result.insert(make_pair(cellIndex, Tree::trees[cellIndex].getSeedPlant()));
            // }
            // cerr << Tree::numberOfTrees << "\n";
        }
        set<pair<int, pair<int, int>>, decltype(cmp2)> result(Tree::getSeedPlant());
    
        for (pair<int, pair<int, int>> r : result) {
            cerr << r.first << " " << r.second.first << " " << r.second.second << endl;
        }
        // for (int i = 0; i < Tree::numberOfTrees; ++i) {
        //     trees[i]
        // }
        int numberOfPossibleMoves;
        cin >> numberOfPossibleMoves; cin.ignore();
        // DEBUGVAR(numberOfPossibleMoves);
        for (int i = 0; i < numberOfPossibleMoves; i++) {
            string possibleMove;
            getline(cin, possibleMove);
        }

        bool tryNext = true;
        if (Tree::countMyTrees(3) >= CMLT_EXP) {
            tryNext = completeTrees();
        }
        if (tryNext && Tree::countMyTrees(0) < 2 && Tree::countMyTrees(-1) < (MAX_TREE * (1 - pow(day / 23.0, 17))) && result.size() > 0) {
            int i = 0;
            tryNext = false;
            cerr << "S" << endl;
            for (const pair<int, pair<int, int>> &item : result) {
                Tree &tree = Tree::trees[item.first];
                Cell &cell = Cell::cells[item.second.first];
                if (cellDistance(*tree.cell, cell) > tree.size) {
                    if (action::grow(tree)) {
                        break;
                    }
                } else {
                    if (action::seed(tree, cell)) {
                        break;
                    }
                }
                ++i;
            }
            if (i == result.size()) {
                // action::wait();
                tryNext = true;
            }
        }
        if (tryNext) {
            cerr << "G" << endl;
            set<pair<int, int>, decltype(cmp)> seeds = Tree::getGrowSeed();
            int i = 0;
            for (const pair<int, int> &growSeed : seeds) {
                if (growSeed.first != -1 && action::grow(Tree::trees[growSeed.first])) {
                    break;
                }
                ++i;
            }
            if (i == seeds.size()) {
                action::wait();
            }
        }


        cerr << "Time difference = " << std::chrono::duration_cast<std::chrono::milliseconds>(NOW - start).count() << "[ms]" << std::endl;

        // GROW cellIdx | SEED sourceIdx targetIdx | COMPLETE cellIdx | WAIT <message>
        // cout << "WAIT" << endl;
    }
}