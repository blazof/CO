#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <unordered_set>
#include <string>
#include <sstream>
#include <algorithm>
#include <climits>
#include <fstream>
#include <stdarg.h>
#include <strings.h>

#include "levenstein.h"

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
        if (i > n - k - 3) {
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

vector<string> positiveErrorHandler(const vector<string>& spectrum, const vector<string>& positiveErrors) {
    vector<string> combinedVector = spectrum;
    for (const auto& positiveError : positiveErrors) {
        combinedVector.push_back(positiveError);
    }
    return combinedVector;
}
void pathByOne(vector<int> &notVisited, vector<string> &spectrum, vector<vector<int>> &graph, int &index, string &output) {
    bool progress = true; // Ensure progress to avoid infinite loops
    while (progress && !notVisited.empty()) {
        progress = false;
        vector<int> pickOnePath;

        // Iterate through the graph's adjacency matrix for the current index
        for (int i = 0; i < spectrum.size(); i++) {

            if (graph[index][i] == 1 && contains(notVisited, i)) {
                pickOnePath.push_back(i);
            }


            if(i==spectrum.size()-1 && !pickOnePath.empty()) {
                int randomPath = rand() % pickOnePath.size();
                string oligo = spectrum[pickOnePath[randomPath]];
                int shorter = min(spectrum[index].size(), spectrum[pickOnePath[randomPath]].size());

                // Append overlapping part to the output
                oligo = oligo.substr(shorter - 1, oligo.size());
                output += oligo;

                // Update `index` and mark vertex as visited
                index = pickOnePath[randomPath];
                //notVisited.erase(find(notVisited.begin(), notVisited.end(), pickOnePath[randomPath]));
                notVisited.erase(find(notVisited.begin(), notVisited.end(), index));
//  tutaj zmiana bo musimy ten index usunac


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
        cout << "2. Algorytm naiwny" << endl;
        cout << "3. Metaheurystyka" << endl;
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
            case 2:
                cout << "Naiwny in progress" << endl;
                break;
            case 3:
                cout << "Metaheurystyka in progress" << endl;
                break;
            default:
                cout << "Żadna z opcji nie jest prawidłowa. Wybierz jeszcze raz." << endl;
                repeat = true;
        }
    } while (repeat);
}
int minDistance(int dist[], bool sptSet[], int V) {
    int min = INT_MAX, min_index;
    for (int v = 0; v < V; v++) {
        if (sptSet[v] == false && dist[v] <= min) {
            min = dist[v], min_index = v;
        }
    }
    return min_index;
}

void dijkstra(const vector<vector<int>> &graph, int src, int V, vector<int> &dist) {
    bool sptSet[V];  // Tablica sprawdzająca, które wierzchołki zostały odwiedzone
    fill(dist.begin(), dist.end(), INT_MAX);  // Inicjalizujemy odległości jako nieskończoność
    fill(sptSet, sptSet + V, false);  // Wszystkie wierzchołki są na początku nieodwiedzone

    dist[src] = 0;  // Odległość do samego siebie wynosi 0

    // Pętla Dijkstry
    for (int count = 0; count < V - 1; count++) {
        int u = minDistance(dist.data(), sptSet, V);  // Wybieramy wierzchołek o najmniejszej odległości
        sptSet[u] = true;  // Oznaczamy go jako odwiedzony

        // Aktualizujemy odległości sąsiadów
        for (int v = 0; v < V; v++) {
            if (!sptSet[v] && graph[u][v] != 0 && dist[u] != INT_MAX && dist[u] + graph[u][v] < dist[v]) {
                dist[v] = dist[u] + graph[u][v];
            }
        }
    }
}

int main() {
    srand(static_cast<unsigned>(time(0)));

    int n = 400, k = 8, delta_k = 2, nError = 0, pError = 0, probablePositive = 0;
    string input;
    bool repAllowed = true;
    string DNA, primer, DNADWA;
    vector<string> idealSpectrum, spectrum, positiveErrors;

    menu(DNA, n, k, delta_k, repAllowed, nError, pError, probablePositive);

    idealSpectrum = generateIdealSpectrum(k, n, DNA, delta_k, repAllowed);
    primer = idealSpectrum[0];

    cout << "Wygenerowane DNA: " << DNA << endl;

    /* DNADWA = generateDNA(n=400, k=8, delta_k=2, repAllowed=true, nError=0, pError=0, probablePositive=0);
    cout<< "Drugie DNA:" << DNADWA << endl;                                                   // TUTAJ JEST TEST MIARY LEVENSTEINA - DZIALA!!!
    cout << "MIARA" << levenshteinDist(DNA, DNADWA) << endl; */

    cout << "Pierwszy oligonukleotyd: " << primer << endl;

    for (const string &element: idealSpectrum) {
        cout << element << " ";
    }
    cout << endl;

    cout << "UWAGA SPEKTRUM!" << endl;
    spectrum = negativeErrorsHandler(idealSpectrum, nError, primer, n, k, delta_k, repAllowed, pError,
                                     probablePositive);
    positiveErrors = positiveErrorGenerator(pError, k, spectrum, delta_k, probablePositive);

    for (const string &element: spectrum) {
        cout << element << " ";
    }
    cout << endl;

    cout << "Tylko pozytywne errory" << endl;
    for (const string &element: positiveErrors) {
        cout << element << " ";
    }
    cout << endl;

    spectrum = positiveErrorHandler(spectrum, positiveErrors);

    cout << "UWAGA SPEKTRUM Z POZYTYWNYMI BŁĘDAMI!" << endl;
    for (const string &element: spectrum) {
        cout << element << " ";
    }
    cout << endl;

    sort(spectrum.begin(), spectrum.end());

    cout << "Posortowane spektrum ze wszystkimi błędami" << endl;
    for (const string &element: spectrum) {
        cout << element << " ";
    }



    //GRAF MOMENT
    vector<vector<int>> graph(spectrum.size(), vector<int>(spectrum.size(), 0));

    for (int i = 0; i < spectrum.size(); i++) {
        for (int j = 0; j < spectrum.size(); j++) {
            if (i == j) {
                continue;
            }
            int size = 0;
            if (spectrum[i].size() > spectrum[j].size()) {
                size = spectrum[j].size();
            } else {
                size = spectrum[i].size();
            }
            //jedynki
            string tmp1 = spectrum[i].substr(spectrum[i].size() - size + 1, size - 1);
            string tmp2 = spectrum[j].substr(0, size - 1);

            if (tmp1 == tmp2) {
                graph[i][j] = 1;
                cout << tmp1 << " == " << tmp2 << endl;
            } else if (nError > 0 && k - delta_k > 2) {

                tmp1 = tmp1.substr(1, tmp1.size());
                tmp2 = tmp2.substr(0, tmp2.size() - 1);
                if (tmp1 == tmp2) {
                    graph[i][j] = 2;
                    cout << tmp1 << " == " << tmp2 << endl;
                } else if (k - delta_k > 3) {
                    tmp1 = tmp1.substr(1, tmp1.size());
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
    int index = 0;
    string output = "";
    vector<int> notVisited;

    // Initialize `notVisited` with all indices
    for (int i = 0; i < spectrum.size(); i++) {
        notVisited.push_back(i);
    }

    // Find the primer in the spectrum and set as the starting point
    for (int i = 0; i < spectrum.size(); i++) {
        if (spectrum[i] == primer) {
            index = i;
            output += spectrum[i];
            notVisited.erase(find(notVisited.begin(), notVisited.end(), i)); // Safe removal
            break;
        }
    }


//    TRZEBA NAPRAWIC TO ZE LOSUJE NAM INDEX NA KTORYM SIE ZATRZYMALISMY I DO NIEGO PRZECHODZI!!!!!!!


    int toVisitPercent = 0.2;
    vector<int> dist(spectrum.size(), INT_MAX);
    pathByOne(notVisited, spectrum, graph, index, output);
    cout << "Reconstructed sequence: " << output << endl;
    cout << "SIZE" << endl;
    cout << notVisited.size() << " " << spectrum.size() << endl;

    float condition = static_cast<float>(notVisited.size()) / spectrum.size();
    cout << "condition = " << condition << endl;

    if (condition > toVisitPercent) {
        cout << "CHUJ" << endl;
        vector<int> toVisit = verticesToVisit(graph, notVisited, spectrum, index);

        // Uruchamiamy Dijkstrę raz, obliczamy ścieżki z `index`
        cout << "Obliczam odległości od wierzchołka: " << index << endl;
        dijkstra(graph, index, spectrum.size(), dist);

        // Wyświetlamy odległości tylko do wierzchołków w `toVisit`
        cout << "Odległości do wierzchołków wylosowanych z toVisit:" << endl;
        int minDist = INT_MAX;
        int nextVertex = -1;

        for (int vertex : toVisit) {
            if (dist[vertex] == INT_MAX) {
                cout << "Do wierzchołka " << vertex << ": brak połączenia (INT_MAX)" << endl;
            } else {
                cout << "Do wierzchołka " << vertex << ": " << dist[vertex] << endl;

                // Szukamy najbliższego wierzchołka
                if (dist[vertex] < minDist) {
                    minDist = dist[vertex];
                    nextVertex = vertex;
                }
            }
        }

        // Wybieramy wierzchołek o najmniejszej odległości
        if (nextVertex != -1) {
            cout << "Najbliższy wierzchołek do odwiedzenia to: " << nextVertex
                 << " z odległością: " << minDist << endl;
        } else {
            cout << "Nie znaleziono żadnego odpowiedniego wierzchołka do odwiedzenia!" << endl;
        }
    }


    // pathByOne(notVisited,spectrum,graph,index,output);
        //cout<<index<<endl;
        //vector<int> toVisit = verticesToVisit(graph,notVisited,spectrum,index);

        cout << "Reconstructed sequence: " << output << endl;
    return 0;
    }
//POPRAWOIC LOSOWANIE (LOSUJE TEN SAM), NIE DODAJE SCIEKZI Z DIKSRY DTAM DALEJ (WIEMY JAK JEST DLUGA ALE JEJ NIE REALIZUJEMY) JESZCZ :)
