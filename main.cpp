#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <unordered_set>
#include <string>
#include <sstream>
#include <algorithm>
#include <cassert>
#include <climits>
#include <fstream>
#include <stdarg.h>
#include <strings.h>
#include "levenstein.h"
#include <tuple>

using namespace std;

template <typename T>
bool contains(const vector<T>& vec, const T& value) {
    return find(vec.begin(), vec.end(), value) != vec.end();
}

string readDNAFromFile(string& DNA, int &n, int &k, int &delta_k, int &nError, int &pError, int &probablePositive, bool &repAllowed) {
    string filename = "DNA.txt";
    ifstream inFile(filename.c_str());

    if (inFile.is_open()) {
        getline(inFile, DNA);
        inFile >> n;
        inFile >> k;
        inFile >> delta_k;
        inFile >> probablePositive;
        inFile >> nError;
        inFile >> pError;
        inFile >> repAllowed;

        inFile.close();
        return DNA;
    } else {
        cerr << "Nie można otworzyć pliku do odczytu: " << filename << endl;
        return "";
    }
}

void saveToFile(const string DNA, const int n, const int k, const int delta_k, const int nError, const int pError, const int probablePositive,const bool repAllowed) {
    string filename = "DNA.txt";
    ofstream outFile(filename.c_str());
    if(outFile.is_open()) {
        cout<<"PLIK OPEN"<<endl;
        outFile << DNA<< "\n"<< n << "\n" << k << "\n" << delta_k << "\n" << probablePositive << "\n"<<nError << "\n" << pError<< "\n" << repAllowed<<"\n";
        outFile.close();
    }else {
        cerr << "Nie można otworzyć pliku do zapisu. " << filename << endl;
    }
}

string generateDNA(int &n, int &k, int &delta_k, bool &repAllowed, int &nError, int &pError, int &probablePositive) {
    string input;

    // Wczytywanie długości łańcucha z domyślną wartością
    cout << "Podaj długość łańcucha (domyślnie " << n << "): ";
    getline(cin, input);
    getline(cin, input);
    if (!input.empty() && isdigit(input[0])) {
        stringstream(input) >> n;
    }

    // Wczytywanie długości oligonukleotydów z domyślną wartością
    cout << "Podaj długość oligonukleotydów (domyślnie " << k << "): ";
    getline(cin, input);
    if (!input.empty() && isdigit(input[0])) {
        stringstream(input) >> k;
    }

    // Wczytywanie delta_k
    cout << "Podaj możliwą zmiennosć długości oligonukleotydów (delta_k) (domyślnie " << delta_k << "): ";
    getline(cin, input);
    if (!input.empty() && isdigit(input[0])) {
        stringstream(input) >> delta_k;
    }

    // Czy powtórzenia są dozwolone?
    cout << "Czy powtórzenia są dozwolone? T/N (domyślnie T): ";
    getline(cin, input);
    if (!input.empty()) {
        if (input == "N") {
            repAllowed = false;
        }else {
            repAllowed = true;
        }
    }

    // Wczytywanie liczby błędów negatywnych z domyślną wartością
    cout << "Podaj ilość błędów negatywnych (domyślnie " << nError << "): ";
    getline(cin, input);
    if (!input.empty() && isdigit(input[0])) {
        stringstream(input) >> nError;
    }

    // Wczytywanie liczby błędów pozytywnych z domyślną wartością
    cout << "Podaj ilość błędów pozytywnych (domyślnie " << pError << "): ";
    getline(cin, input);
    if (!input.empty() && isdigit(input[0])) {
        stringstream(input) >> pError;
    }

    cout << "Czy błędy pozytywne mają być realistyczne? T/N (domyślnie T): ";
    getline(cin, input);
    if (!input.empty()) {
        if (input == "T") {
            probablePositive =1;
        }else {
            probablePositive =0;
        }
    }

    string DNA;
    for (int i = 0; i < n; i++) {
        const char nucleotides[] = {'A', 'C', 'T', 'G'};
        char generatedNucleotide = nucleotides[rand() % 4];
        DNA += generatedNucleotide;
    }
    //zapis do pliku
    saveToFile(DNA,n,k,delta_k,nError,pError,probablePositive, repAllowed);
    return DNA;
}

vector<string> generateIdealSpectrum(const int k, const int n, string& DNA, const int delta_k, bool repAllowed) {
    int shift = 0;
    vector<string> idealSpectrum;
    string oligonucleotide = "";

    for (int i = 0; i <= n - k; i++) {
        if ((i > n - k - 3) || i==0) {
            oligonucleotide = DNA.substr(i, k);
        } else {
            if (delta_k > 0) {
                shift = rand() % (delta_k + 1);
                if (rand() % 2 == 0) {
                    shift *= -1;
                }
            }

            oligonucleotide = DNA.substr(i, k + shift);


        }

        if(!repAllowed) {
            if(contains(idealSpectrum, oligonucleotide)) {
                cout<<"POWTÓRZENIE"<<endl;
                do {
                    const char nucleotides[] = {'A', 'C', 'T', 'G'};
                    char generatedNucleotide = nucleotides[rand() % 4];
                    int randomIndex = rand() % oligonucleotide.length();
                    cout << "OLIGO PRZED ZMIANA: " << oligonucleotide <<endl;
                    oligonucleotide[randomIndex] = generatedNucleotide;
                    cout<<DNA<<endl;
                    cout << "TU BYŁO POWTÓRZENIE: " << oligonucleotide << " DNA: " << DNA[i+randomIndex]<< endl;

                    DNA[i+randomIndex] = generatedNucleotide;
                    cout<<DNA<<endl;
                }while (contains(idealSpectrum,oligonucleotide));

            }
        }
        idealSpectrum.push_back(oligonucleotide);
    }
    return idealSpectrum;
}

vector<string> negativeErrorsHandler(const vector<string>& spectrum, const int nError, const string& primer, int &n, int &k, int &delta_k, bool &repAllowed, int &pError, int &probablePositive) {
    int repeats = 0, difference = 0;

    // Tworzymy zbiór, który automatycznie usuwa duplikaty
    unordered_set<string> uniqueSet(spectrum.begin(), spectrum.end());

    // Zliczamy powtórzenia
    repeats = spectrum.size() - uniqueSet.size();
    vector<string> uniqueVec(uniqueSet.begin(), uniqueSet.end());

    difference = nError - repeats;
    cout<<"Powtórzenia: "<< repeats <<endl;
    /*
    cout << "Powtórzenia: " << repeats << "Czy chcesz kontynować (w przeciwnym razie instancja będzie generowana na nowo. T/N (domyślnie T)" << endl;

    string input;
    getline(cin, input);

    if (!input.empty()) {
        if (input == "N") {
            // generujemy na nowo

        }
    }else {

    }
    */

    for (int i = 0; i < difference; i++) {
        if (!uniqueVec.empty()) {
            while (true) {
                int position = rand() % uniqueVec.size();
                if (uniqueVec[position] != primer) {
                    uniqueVec.erase(uniqueVec.begin() + position);
                    break;
                }
            }
        }
    }

// ACTCTTTCACGGACTACACACTGTGGCAGCGATTCCAGTACCGCGCCAAGCGCACAGTAGAAGTTAGTACGCGTTTAGTCTGAGATAGAATAACGTTTCTAAGGTCCTAGAGGATGTTTGGCGGTTTACGTCGCACTCAGAATGAATGGCTTTCCAAGCCTGTCCTGGGGTGGAGGCCTATATATGCTGTGCCCGTGTTTCACAAAGTTTTATACAGTTTACACTGACCTGTTCGGCCCGGGATGCTGGCCGTCACAATCTTGCGCAAGTATCAACCGAAACGGGGGGATCTTAGGTAGC
    return uniqueVec;
}

vector<int> verticesToVisit(const vector<vector<int>> &graph, vector<int> &notVisited, vector<string> spectrum, int &finalIndex) {
    vector<int> vertices;       // Lista wierzchołków do odwiedzenia (indexy)
    vector<int> toVisit = notVisited;  // Kopia listy wierzchołków, które jeszcze nie zostały odwiedzone
    int index = 0;

    // Chcemy odwiedzić 3 wierzchołki
    for (int i = 0; i < 10; i++) {
        bool found = false;  // Flaga wskazująca, czy znaleziono wierzchołek do odwiedzenia
        index = rand() % toVisit.size();  // Losujemy losowy indeks wierzchołka

        for (int j = 0; j < spectrum.size(); j++) {
            // Sprawdzamy sąsiadów wylosowanego wierzchołka
            if (graph[index][j] != 0 && contains(toVisit, j) && index != finalIndex) {  // Sprawdzamy, czy sąsiad jest nieodwiedzony
                vertices.push_back(index);  // Dodajemy do listy odwiedzonych
                // Usuwamy wierzchołek z toVisit, używając podejścia z indeksem
                for (int x = 0; x < toVisit.size(); x++) {
                    if (x == index) {
                        toVisit.erase(std::remove(toVisit.begin(), toVisit.end(), x), toVisit.end());                        break;  // Przerywamy po usunięciu
                    }
                }
                found = true;
                break;  // Przerywamy pętlę, gdy znajdziemy sąsiada, którego jeszcze nie odwiedziliśmy
            }
        }

        if (!found) {
            i--;  // Jeśli nie znaleziono wierzchołka do odwiedzenia, powtarzamy iterację
        }
    }

    // Debug: Wypisanie wierzchołków do odwiedzenia
    cout << "Elementy do odwiedzenia: ";
    for (int element : vertices) {
        cout << element << " ";
    }
    cout << endl;

    return vertices;
}

vector<string> positiveErrorGenerator(const int pError, const int k, const vector<string>& spectrum, const int delta_k, const int probablePositive) {
    vector<string> positiveErrors;

    for (int i = 0; i < pError; i++) {
        string positiveError;
        do {
            int shift=0;
            if (delta_k > 0) {
                shift = rand() % (delta_k + 1);
                if (rand() % 2 == 0) {
                    shift *= -1;
                }
            }
            positiveError = "";
            for (int j = 0; j < k+shift; j++) {
                const char nucleotides[4] = {'A', 'C', 'T', 'G'};
                const char generatedNucleotide = nucleotides[rand() % 4];
                positiveError += generatedNucleotide;
            }
        } while (contains(spectrum, positiveError) || contains(positiveErrors, positiveError));
        positiveErrors.push_back(positiveError);
    }
    return positiveErrors;
}

bool hasUnvisitedAdj(vector<vector<int>> updatedGraph, int V, int index, vector<int> &notVisited){
    for(int i =0; i < V;i++){
        if(updatedGraph[index][i] != INT_MAX && find(notVisited.begin(), notVisited.end(), i)!= notVisited.end()){
            return true;
        }
    }

    return false;
}

vector<string> positiveErrorHandler(const vector<string>& spectrum, const vector<string>& positiveErrors) {
    vector<string> combinedVector = spectrum;
    for (const auto& positiveError : positiveErrors) {
        combinedVector.push_back(positiveError);
    }
    return combinedVector;
}

void pathByOne(vector<int> &notVisited, vector<string> &spectrum, vector<vector<int>> &graph, int &index, string &output, vector<vector<int>> &updatedGraph, vector<int> &seq) {
    bool progress = true; // Ensure progress to avoid infinite loops
    while (progress && !notVisited.empty()) {
        progress = false;
        vector<int> pickOnePath;

        // Iterate through the graph's adjacency matrix for the current index
        for (int i = 0; i < spectrum.size(); i++) {
            if (index < graph.size() && i < graph[index].size() && graph[index][i] == 1 &&
                contains(notVisited, i) && hasUnvisitedAdj(updatedGraph, spectrum.size(), i, notVisited)) {
                pickOnePath.push_back(i);
            }

            if (i == spectrum.size() - 1 && !pickOnePath.empty()) {
                int randomPath = rand() % pickOnePath.size();

                // Safeguard: Check if pickOnePath is non-empty
                assert(!pickOnePath.empty());

                // Append overlapping part to the output
                string oligo = spectrum[pickOnePath[randomPath]];
                int shorter = min(spectrum[index].size(), spectrum[pickOnePath[randomPath]].size());

                // Safeguard: Ensure shorter > 0 before using substr
                if (shorter > 0) {
                    oligo = oligo.substr(shorter - 1, oligo.size());
                    output += oligo;
                }

                // Update index and mark vertex as visited
                index = pickOnePath[randomPath];

                // Safeguard: Ensure index exists in notVisited before erasing
                auto it = find(notVisited.begin(), notVisited.end(), index);
                if (it != notVisited.end()) {
                    notVisited.erase(it);
                }

                seq.push_back(index);

                // Update the adjacency matrix
                for (int x = 0; x < spectrum.size(); x++) {
                    if (index < updatedGraph.size() && x < updatedGraph[index].size()) {
                        updatedGraph[index][x] = 0;
                        updatedGraph[x][index] = 0;
                    }
                }

                progress = true; // Continue to the next vertex
                pickOnePath.clear();
                break; // Exit the inner loop to process the next vertex
            }
        }
    }
}

void menu(string &DNA, int &n, int &k, int &delta_k, bool &repAllowed, int &nError, int &pError, int &probablePositive) {
    bool repeat = false;
    do {
        int choice = 0;
        repeat = false;

        cout << "     Menu główne" << endl;
        cout << "1. Generator instancji" << endl;
        cin >> choice;

        switch (choice) {
            case 1:
                cout << "1. Wczytaj DNA z pliku" << endl;
                cout << "2. Generuj ręcznie" << endl;
                cin >> choice;

                switch (choice) {
                    case 1:
                        DNA = readDNAFromFile(DNA,n, k, delta_k, nError, pError, probablePositive, repAllowed);
                        break;
                    case 2:
                        DNA = generateDNA(n = 400, k = 8, delta_k = 2, repAllowed = true, nError = 0, pError = 0, probablePositive = 0);
                        break;
                    default:
                        cout << "Podałeś złą opcję menu, wybierz jeszcze raz." << endl;
                        repeat = true;
                }
                break;
            default:
                cout << "Żadna z opcji nie jest prawidłowa. Wybierz jeszcze raz." << endl;
                repeat = true;
        }
    } while (repeat);
}

using namespace std;

vector<vector<int>> generateGraph(const vector<string>& spectrum, const int delta_k, const int k) {
    vector<vector<int>> graph(spectrum.size(), vector<int>(spectrum.size(), 0));

    for (int i = 0; i < spectrum.size(); i++) {
        for (int j = 0; j < spectrum.size(); j++) {
            if (i == j) {
                continue;
            }
            int size = min(spectrum[i].size(), spectrum[j].size());
            string tmp1 = spectrum[i].substr(spectrum[i].size() - size + 1, size - 1);
            string tmp2 = spectrum[j].substr(0, size - 1);

            if (tmp1 == tmp2) {
                graph[i][j] = 1;
                cout << tmp1 << " == " << tmp2 << endl;
            } else if (k - delta_k > 2) {
                tmp1 = tmp1.substr(1, tmp1.size() - 1);
                tmp2 = tmp2.substr(0, tmp2.size() - 1);
                if (tmp1 == tmp2) {
                    graph[i][j] = 2;
                    cout << tmp1 << " == " << tmp2 << endl;
                } else if (k - delta_k > 3) {
                    tmp1 = tmp1.substr(1, tmp1.size() - 1);
                    tmp2 = tmp2.substr(0, tmp2.size() - 1);
                    if (tmp1 == tmp2) {
                        graph[i][j] = 3;
                    }
                }
            }
        }
    }

    for (int i = 0; i < spectrum.size(); i++) {
        for (int j = 0; j < spectrum.size(); j++) {
            cout << graph[i][j] << " ";
        }
        cout << endl;
    }
    return graph;
}

string startPath(const vector<string>& spectrum, vector<int>& notVisited, const string primer, int& index, vector<vector<int>>& updatedGraph) {
    string output = "";
    for (int i = 0; i < spectrum.size(); i++) {
        notVisited.push_back(i);
    }

    for (int i = 0; i < spectrum.size(); i++) {
        if (spectrum[i] == primer) {
            index = i;
            output += spectrum[i];
            notVisited.erase(find(notVisited.begin(), notVisited.end(), i)); // Safe removal

            for (int x = 0; x < spectrum.size(); x++) {
                updatedGraph[index][x] = 0;
                updatedGraph[x][index] = 0;
            }
            break;
        }
    }
    return output;
}

vector<int> distInit(vector<vector<int>> updatedGraph, int V, int index, vector<int>& parent) {
    vector<int> dist(V, INT_MAX);
    for (int i = 0; i < V; i++) {
        parent[i] = i;
    }
    dist[index] = 0;
    return dist;
}

int getNearest(int V, vector<int> dist, vector<bool> visited) {
    int minValue = INT_MAX;
    int minNode = -1;
    for (int i = 0; i < V; i++) {
        if (!visited[i] && dist[i] < minValue) {
            minValue = dist[i];
            minNode = i;
        }
    }
    return minNode;
}

void dijkstra(vector<vector<int>> updatedGraph, int V, vector<int>& dist, vector<bool>& visited, vector<int>& parent) {
    for (int i = 0; i < V; i++) {
        int nearest = getNearest(V, dist, visited);
        if (nearest == -1) break; // If no more reachable nodes
        visited[nearest] = true;

        for (int adj = 0; adj < V; adj++) {
            if (updatedGraph[nearest][adj] != INT_MAX && dist[adj] > dist[nearest] + updatedGraph[nearest][adj]) {
                dist[adj] = dist[nearest] + updatedGraph[nearest][adj];
                parent[adj] = nearest;
            }
        }
    }
}

void displayDistances(vector<int>& dist, vector<int>& parent, int V, int start) {
    cout << "Źródło: " << start << endl;
    cout << "Wierzchołek\tOdległość\tŚcieżka" << endl;

    for (int i = 0; i < V; i++) {
        if (dist[i] == INT_MAX) {
            cout << i << "\t\tNIESKOŃCZONE\t-" << endl;
        } else {
            cout << i << "\t\t" << dist[i] << "\t\t" << i;
            int p = parent[i];
            while (p != start && p != parent[p]) { // Checking to avoid looping back
                cout << " <- " << p;
                p = parent[p];
            }
            if (p == start) {
                cout << " <- " << start;
            }
            cout << endl;
        }
    }
}

string mergeSequences(const std::string& seq1, const std::string& seq2) {
    int maxOverlap = 0;
    int n1 = seq1.size();
    int n2 = seq2.size();

    for (int i = 1; i <= std::min(n1, n2); i++) {
        if (seq1.substr(n1 - i) == seq2.substr(0, i)) {
            maxOverlap = i;
        }
    }

    return seq1 + seq2.substr(maxOverlap);
}

vector<vector<float>> matrixACO(int V, vector<vector<int>> rankedVertices, vector<float> values, int rankMatrix) {
    vector<vector<float>> matrix(V, vector<float>(V, 0));
    for (int i = 0; i < rankMatrix; i++) {
        for (int j = 0; j < rankedVertices[i].size() - 1; j++) {
            if (rankedVertices[i][j] < V && rankedVertices[i][j + 1] < V) {
                matrix[rankedVertices[i][j]][rankedVertices[i][j + 1]] += values[i];
            } else {
                cerr << "Index out of bounds in matrixACO: " << rankedVertices[i][j] << " or " << rankedVertices[i][j + 1] << endl;
            }
        }
    }
    for (int i = 0; i < V; i++) {
        for (int j = 0; j < V; j++) {
            cout << matrix[i][j] << " ";
        }
        cout << endl;
    }
    return matrix;
}

vector<int> greedyAlgorithm(vector<vector<int>> updatedGraph, int V, vector<vector<int>> graph, vector<string> &spectrum, int index, vector<int> notVisited, string output, string DNA, vector<int> seq) {
    vector<int> sequence = seq;
    pathByOne(notVisited, spectrum, graph, index, output, updatedGraph,sequence);

    while (static_cast<float>(notVisited.size()) / spectrum.size() > 0.2) {
        vector<bool> visited(V,false);
        vector<int> parent(V, 0);
        vector<int> dist(V, 0);

        for (int x = 0; x < spectrum.size(); x++) {
            updatedGraph[index][x] = graph[index][x];
            updatedGraph[x][index] = graph[x][index];
        }

        for (int i = 0; i < spectrum.size(); i++) {
            for (int j = 0; j < spectrum.size(); j++) {
                if (updatedGraph[i][j] == 0) {
                    updatedGraph[i][j] = INT_MAX;
                }
            }
        }

        cout << endl;


        dist = distInit(updatedGraph, V, index, parent);
        dijkstra(updatedGraph, V, dist, visited, parent);

        /*   cout << "Final distances from source:" << endl;
           for (int i = 0; i < dist.size(); i++) {
               if (dist[i] == INT_MAX) cout << i << ": INF" << endl;
               else cout << i << ": " << dist[i] << " - parent: " << parent[i] << endl;
           }

           displayDistances(dist, parent, spectrum.size(), index);
   */
        if (notVisited.empty()) break;

        int minValue = INT_MAX;
        int vertex = -1;
        vector<int> toMerge;
        int counter=0;
        bool end=false;
        for (int i = 0; i < 50; i++) {
            counter++;
            if(counter > spectrum.size()){
                cout<<"KONIEC"<<endl;
                end = true;
                break;
            }
            int randomIndex = rand() % spectrum.size();
            if (randomIndex == index || dist[randomIndex] == INT_MAX || !hasUnvisitedAdj(updatedGraph,V, randomIndex, notVisited)) {
                i--;
                continue;
            }
            if (dist[randomIndex] < minValue) {
                minValue = dist[randomIndex];
                vertex = randomIndex;
            }
        }
        if(end){
            break;
        }
        cout << "Selected vertex: " << vertex << endl;
        toMerge.push_back(vertex);

        if (minValue == INT_MAX) {
            cout << "No reachable vertex found" << endl;
        } else {
            int p = parent[vertex];
            while (p != index && p != parent[p]) {
                cout << " <- " << p;
                toMerge.push_back(p);
                p = parent[p];
            }
            if (p == index) {
                cout << " <- " << index;
            }
        }
        if(toMerge.size() > 0) {
            reverse(toMerge.begin(), toMerge.end());

            for (int i = 0; i < toMerge.size(); i++) {
                output = mergeSequences(output, spectrum[toMerge[i]]);
                auto it = std::find(notVisited.begin(), notVisited.end(), toMerge[i]);
                if (it != notVisited.end()) {
                    notVisited.erase(it);
                }
                sequence.push_back(toMerge[i]);
            }
            index = toMerge[0];

        }

        cout << endl;

        cout << "Output after Dijkstra: " << output << endl;


        if(notVisited.empty()) {
            break;
        }

        pathByOne(notVisited, spectrum, graph, index, output, updatedGraph,sequence);
        cout << "Output after path by one: " << output << endl;

    }
    cout<<endl;

    if (output.length() > DNA.length()) {
        output = output.substr(0, DNA.length());
    }

    cout<<"DNA x: "<<DNA<<endl;
    cout<<"Final: "<<output <<endl;

    cout<<"Lev: "<<levenshteinDist(DNA,output)<<endl;
    return sequence;
}

string rankingACO(vector<vector<float>> &finisedMatrix,int param, int rankMatrix, int V, string output, vector<vector<int>> updatedGraph, vector<vector<int>> graph, vector<string> &spectrum, int index, vector<int> notVisited, string &DNA, vector<int> seq) {
    vector<float> values = {1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1};
    vector<vector<float>> matrix(V, vector<float>(V, 0));
    vector<vector<int>> rankedVertices;
    vector<int> ones(param, 0); // param - ile wierzcholkow jest w rankingu

    for (int i = 0; i < param; i++) {
        rankedVertices.push_back(greedyAlgorithm(updatedGraph, V, graph, spectrum, index, notVisited, output, DNA, seq));
        for (int j = 0; j < rankedVertices[i].size() - 1; j++) {
            if (rankedVertices[i][j] < V && rankedVertices[i][j + 1] < V) {
                if (graph[rankedVertices[i][j]][rankedVertices[i][j + 1]] == 1) {
                    ones[i]++;
                }
            } else {
                cerr << "Index out of bounds in ACO: " << rankedVertices[i][j] << " or " << rankedVertices[i][j + 1] << endl;
            }
        }
    }

    // Tworzenie wektora par (ones, {rankedVertices, size})
    vector<tuple<int, vector<int>, size_t>> paired;
    for (int i = 0; i < param; i++) {
        paired.push_back(make_tuple(ones[i], rankedVertices[i], rankedVertices[i].size()));
    }

    // Sortowanie par według wartości `ones` malejąco, a następnie według rozmiaru `rankedVertices`
    sort(paired.begin(), paired.end(), [](const tuple<int, vector<int>, size_t> &a, const tuple<int, vector<int>, size_t> &b) {
        if (get<0>(a) == get<0>(b)) {
            return get<2>(a) > get<2>(b);
        }
        return get<0>(a) > get<0>(b);
    });

    // Rozdzielenie posortowanych par z powrotem na `rankedVertices` i `ones`
    for (int i = 0; i < param; i++) {
        ones[i] = get<0>(paired[i]);
        rankedVertices[i] = get<1>(paired[i]);
    }

    // Wyświetlanie posortowanych wyników
    for (int i = 0; i < param; i++) {
        for (int j = 0; j < rankedVertices[i].size(); j++) {
            cout << rankedVertices[i][j] << " ";
        }
        cout << endl;
        cout << ones[i] << endl;
        cout << rankedVertices[i].size() << endl;
    }

    matrix = matrixACO(V, rankedVertices, values, rankMatrix);
    finisedMatrix = matrix;
    return output;
}

//jako index nalezy przekazac primer!!!!
void initACO(string DNA,int ants, int smoothing, int interations, float firstDrawPercentage, int index, vector<vector<float>> matrix, vector<vector<int>> graph, int V, vector<string> spectrum, int n) {
    float drawPercentage = firstDrawPercentage;
    vector<vector<int>> updatedGraph = graph;

    // Rozwiązania dla wszystkich mrówek
    vector<vector<int>> solutions(ants); // Tworzymy wektor z 'ants' pustymi wektorami wewnętrznymi
    vector<string> outputs;
    // Iteracje algorytmu
    for (int j = 0; j < interations; j++) {
        drawPercentage = firstDrawPercentage; // Resetowanie wartości drawPercentage dla każdej iteracji
        // Reset grafu dla nowej iteracji

        // Przejście każdej mrówki
        for (int i = 0; i < ants; i++) {
            string output = spectrum[index];
            bool con = true; // Flaga warunku zakończenia dla mrówki
            int currentIndex = index; // Kopia początkowego wierzchołka dla każdej mrówki
            updatedGraph = graph;
            solutions[i].clear(); // Czyścimy rozwiązanie mrówki (dla nowej dużej iteracji)
            solutions[i].push_back(currentIndex); // Dodanie początkowego wierzchołka

            while (con) {
                int drawValue = rand() % 100;

                if (drawValue > drawPercentage) {
                    // Znalezienie wszystkich możliwych ścieżek z bieżącego wierzchołka
                    vector<int> possiblePaths;

                    for (int x = 0; x < updatedGraph[currentIndex].size(); x++) {
                        if (updatedGraph[currentIndex][x] == 1) { // Szukamy ścieżek o wartości 1
                            for (int y = 0; y < updatedGraph[x].size(); y++) {
                                if (updatedGraph[x][y] != 0) {
                                    possiblePaths.push_back(x);
                                    break;
                                }
                            }
                        }
                    }

                    // Jeśli nie ma ścieżek o wartości 1, szukamy innych ścieżek
                    if (possiblePaths.empty()) {
                        for (int x = 0; x < updatedGraph[currentIndex].size(); x++) {
                            if (updatedGraph[currentIndex][x] != 0) {
                                for (int y = 0; y < updatedGraph[x].size(); y++) {
                                    if (updatedGraph[x][y] != 0) {
                                        possiblePaths.push_back(x);

                                        break;
                                    }
                                }
                            }
                        }
                    }

                    // Jeśli nadal nie ma dostępnych ścieżek
                    if (possiblePaths.empty()) {
                        cout << "Brak ścieżek, tworzymy nowe połączenie!" << endl;

                        bool repeat = true;
                        while (repeat) {
                            repeat = false;
                            int randomVertex = rand() % V;

                            if (contains(solutions[i], randomVertex)) {
                                repeat = true; // Jeśli już odwiedziliśmy ten wierzchołek, losujemy dalej
                            } else {
                                updatedGraph[currentIndex][randomVertex] = 1; // Tworzymy sztuczne połączenie
                                updatedGraph[randomVertex][currentIndex] = 1; // Graf nieskierowany
                                currentIndex = randomVertex;
                                solutions[i].push_back(currentIndex);
                                output = mergeSequences(output, spectrum[currentIndex]);
                                cout << "Nowe połączenie: " << currentIndex << endl;
                            }
                        }

                    } else {
                        // Losowy wybór ścieżki z dostępnych
                        int random = rand() % possiblePaths.size();
                        int randomIndex = possiblePaths[random];
                        // Usuwamy krawędzie, aby nie wracać

                        for(int x=0; x < updatedGraph[currentIndex].size(); x++) {
                            updatedGraph[currentIndex][x] = 0;
                        }
                        currentIndex = randomIndex; // Przechodzimy do nowego wierzchołka
                        solutions[i].push_back(currentIndex); // Dodajemy nowy wierzchołek do rozwiązania
                        output = mergeSequences(output, spectrum[currentIndex]);
                        if(output.size() >= n) {
                        con =false;
                        }
                    }
                }else {
                    vector<float> rouletteValues;
                    vector<int> vertices;

                    for(int x = 0; x < matrix[currentIndex].size(); x++) {
                        if(matrix[currentIndex][x] != 0 &&  updatedGraph[currentIndex][x] != 0) {
                            rouletteValues.push_back(matrix[currentIndex][x]*10);
                            vertices.push_back(x);
                        }
                    }

                    if(!rouletteValues.empty()) {
                    vector<int> sum;
                    sum.push_back(0);

                    cout<<"indexy tych wierzchilkkow"<<endl;
                    for(int x=0; x<vertices.size(); x++) {
                        cout<<vertices[x]<<" ";
                    }cout<<endl;

                    cout<<"warotsci stworzenia do ruletki"<<endl;
                    for(int x=0; x<rouletteValues.size(); x++) {
                        cout<<rouletteValues[x]<<" ";
                    }cout<<endl;



                    cout<< "tu byla macierz wybrana"<<endl;

                    //wygladzanie wartosci

                        // szukanie najwiekszego:
                        float biggestIndex = -1;
                        for(int x=0; x<rouletteValues.size(); x++) {
                            if(rouletteValues[x] > rouletteValues[biggestIndex]) {
                                biggestIndex = x;
                            }
                        }

                            for(int x =0; x<rouletteValues.size(); x++) {
                                float difference = rouletteValues[biggestIndex]-rouletteValues[x];
                                if(difference > smoothing) {
                                    rouletteValues[x] *=2;
                                    rouletteValues[biggestIndex] *= 0.9;
                                }


                            }
                        cout<<"warotsci stworzenia do ruletki"<<endl;
                        for(int x=0; x<rouletteValues.size(); x++) {
                            cout<<rouletteValues[x]<<" ";
                        }cout<<endl;


                        sum.push_back(rouletteValues[0]);
                        for(int x=1; x<vertices.size(); x++) {
                            sum.push_back(rouletteValues[x] + sum[x]);
                        }


                        cout<<"Sumy ruletka done"<<endl;
                        for(int x=0; x<sum.size(); x++) {
                            cout<<sum[x]<<" ";
                        }cout<<endl;

                    //tutaj szukamy wartosci ktora bedzie miedzy przedzialami bierzemy ostatni index
                    int rouletteValue = 1 + rand() % sum[sum.size()-1];

                    for(int x=1; x<sum.size(); x++) {
                        if(rouletteValue <= sum[x] && rouletteValue > sum[x-1]) {
                        // tutaj nowym indeksem bedzie vertices[x]
                            int nextIndex = vertices[x-1];
                            cout<<"Roullete value: "<<rouletteValue<<endl;
                            cout<<"nast index"<<" "<<endl;
                            cout<<nextIndex<<" "<<endl;
                            solutions[i].push_back(nextIndex);
                            currentIndex = nextIndex;
                           mergeSequences(output, spectrum[nextIndex]);
                            for(int x=0; x < updatedGraph[currentIndex].size(); x++) {
                                updatedGraph[currentIndex][x] = 0;
                            }
                        }
                    }
                    }else {
                        cout<<"Brak wierzcholkow z macierzy :<<"<<endl;
                    }
                }

                //tutaj ile jest porytego grafu - potrzebne do zwiekszania prawopodobienstwa wyboru macierzy
                float percentCovered = output.size() / n;
                drawPercentage = percentCovered*1.5;
            }
            outputs.push_back(output);
            // Wypisanie aktualnej ścieżki mrówki
            for (int v = 0; v < solutions[i].size(); v++) {
                cout << solutions[i][v] << " ";
            }
            cout << endl;
            cout<< outputs[i] << endl;
            cout<<"Lev"<<levenshteinDist(DNA,outputs[i]);

            cout << endl;


        }
    }
}


void secondMenu(vector<vector<int>> updatedGraph, int &V, vector<vector<int>> graph, vector<string> &spectrum, int index, vector<int> notVisited, string output, string &DNA, vector<int> &seq, int n) {

    vector<vector<int>> vertexSeq;
    while (1){

        vector<int> sequence = seq;
        cout<<"WYgenerowałeś instancję? ŚWIETNIE. Tu masz 2 menu: "<<endl;
        cout<<"1. Algorytm naiwny"<<endl;
        cout<<"2. Metaheurystyka"<<endl;
        cout<<"3. Wyłącz"<<endl;

        int choice;
        bool repeat = false;
        int ones=0;

        cin>>choice;
        switch(choice) {
            case 1: {
                sequence = greedyAlgorithm(updatedGraph,V,graph,spectrum,index,notVisited,output, DNA,sequence);
                for(int i=0; i<sequence.size(); i++) {
                    cout<<sequence[i]<<" ";
                }cout<<endl;

                for(int i=0; i<sequence.size()-1; i++) {
                    if(graph[sequence[i]][sequence[i+1]]==1) {
                        ones++;
                    }
                }
                cout<<"Ilość 1-ek: "<< ones<<endl;
                cout<<"wielkosc seq" << sequence.size()<<endl;;
                break;
            }
            case 2: {
                vector<vector<float>> matrix(V, vector<float>(V, 0));
                cout<< rankingACO(matrix,20,10,V,output,updatedGraph,graph,spectrum,index,notVisited, DNA,sequence)<<endl;
                initACO(DNA,10,20,1,10,index,matrix,graph, V,spectrum, n);
                break;
            }
            case 3: {
                return;
                break;
            }
            default: {
                cout<<"Zły wybór"<<endl;
            }
        }

    }
}

int main() {
    srand(static_cast<unsigned>(time(0)));

    int n = 400, k = 8, delta_k = 2, nError = 0, pError = 0, probablePositive = 0;
    string input;
    bool repAllowed = true;
    string DNA, primer, DNADWA, output;
    vector<string> idealSpectrum, spectrum, positiveErrors;
    vector<vector<int>> graph, updatedGraph;
    vector<int> notVisited, seq;
    int index = 0;

    menu(DNA, n, k, delta_k, repAllowed, nError, pError, probablePositive);
    idealSpectrum = generateIdealSpectrum(k, n, DNA, delta_k, repAllowed);
    primer = idealSpectrum[0];

    spectrum = negativeErrorsHandler(idealSpectrum, nError, primer, n, k, delta_k, repAllowed, pError, probablePositive);
    positiveErrors = positiveErrorGenerator(pError, k, spectrum, delta_k, probablePositive);
    spectrum = positiveErrorHandler(spectrum, positiveErrors);

    cout << "DNA: " << DNA << endl;
    cout << "primer: " << primer << endl;
    cout << "liczba elementow spektrum: " << spectrum.size() << endl;
    sort(spectrum.begin(), spectrum.end());

    int V = spectrum.size();
    for (int i = 0; i < spectrum.size(); i++) {
        cout << spectrum[i] << endl;
        if(spectrum[i]==primer) {
            seq.push_back(i);
        }
    }

    graph = generateGraph(spectrum, delta_k, k);
    updatedGraph = graph;
    output = startPath(spectrum, notVisited, primer, index, updatedGraph);

    cout << "Reconstructed sequence: " << output << endl;

    for (auto e : notVisited) {
        cout << e << " ";
    }
    cout << endl;
    cout << "z grafu" << endl;

    //vertexSequences.push_back( greedyAlgorithm(updatedGraph,V,graph,spectrum,index,notVisited,output, DNA,seq));
    secondMenu(updatedGraph,V,graph,spectrum,index,notVisited,output, DNA,seq,n);

    return 0;
}
//Nowy comment