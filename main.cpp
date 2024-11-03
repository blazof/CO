#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <unordered_set>
#include <string>
#include <sstream>
#include <algorithm>
#include "levenstein.h"

using namespace std;

template <typename T>
bool contains(const vector<T>& vec, const T& value) {
    return find(vec.begin(), vec.end(), value) != vec.end();
}

string readDNAFromFile(int &n, int &k, int &delta_k, bool &repAllowed, int &nError, int &pError, int &probablePositive) {
    string DNA = "";
    return DNA;
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
        repAllowed = (input == "T");
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

    string DNA;
    for (int i = 0; i < n; i++) {
        const char nucleotides[] = {'A', 'C', 'T', 'G'};
        char generatedNucleotide = nucleotides[rand() % 4];
        DNA += generatedNucleotide;
    }
    return DNA;
}

vector<string> generateIdealSpectrum(const int k, const int n, const string& DNA, const int delta_k) {
    int shift =0;
    vector<string> idealSpectrum;
    string oligonucleotide = "";




    for (int i = 0; i <= n - k; i++) {

        if(i > n-k -3) {
            oligonucleotide = DNA.substr(i, k);
        }else {
            if (delta_k > 0) {
                shift = rand() % (delta_k+1);
                if (rand() % 2 == 0) {
                    shift *= -1;
                }
            }
             oligonucleotide = DNA.substr(i, k+shift);
        }

        idealSpectrum.push_back(oligonucleotide);
    }
    return idealSpectrum;
}

vector<string> negativeErrorsHandler(const vector<string>& spectrum, const int nError, const string& primer) {
    int repeats = 0, difference = 0;

    // Tworzymy zbiór, który automatycznie usuwa duplikaty
    unordered_set<string> uniqueSet(spectrum.begin(), spectrum.end());

    // Zliczamy powtórzenia
    repeats = spectrum.size() - uniqueSet.size();
    vector<string> uniqueVec(uniqueSet.begin(), uniqueSet.end());

    difference = nError - repeats;

    cout << "Powtórzenia: " << repeats << endl;

    if (difference > 0) {
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
    } else if (difference < 0) {
        cout << "Nie obsługujemy tego jeszcze :3" << endl;
    }

    return uniqueVec;
}

vector<string> positiveErrorGenerator(const int pError, const int k, const vector<string>& spectrum, const int delta_k, const int probablePositive) {
    vector<string> positiveErrors;

    for (int i = 0; i < pError; i++) {
        string positiveError;

            do {
                positiveError = "";
                for (int j = 0; j < k ; j++) {
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
                        DNA = readDNAFromFile(n, k, delta_k, repAllowed, nError, pError, probablePositive);
                        break;
                    case 2:
                        DNA = generateDNA(n, k, delta_k, repAllowed, nError, pError, probablePositive);
                        cout << "2. Ręcznie" << endl;
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

int main() {
    srand(static_cast<unsigned>(time(0)));
    int n = 400, k = 8, delta_k = 2, nError = 0, pError = 0, probablePositive = 0;
    string input;
    bool repAllowed = true;
    string DNA, primer, DNADWA;
    vector<string> idealSpectrum, spectrum, positiveErrors;

    menu(DNA, n, k, delta_k, repAllowed, nError, pError, probablePositive);

    idealSpectrum = generateIdealSpectrum(k, n, DNA, delta_k);
    primer = idealSpectrum[0];

    cout << "Wygenerowane DNA: " << DNA << endl;

    /* DNADWA = generateDNA(n, k, delta_k, repAllowed, nError, pError, probablePositive);
    cout<< "Drugie DNA:" << DNADWA<<endl;                                                     TUTAJ JEST TEST MIARY LEVENSTEINA - DZIALA!!!
    cout<< "MIARA" << levenshteinDist(DNA, DNADWA)<<endl; */

    cout << "Pierwszy oligonukleotyd: " << primer << endl;

    for (const string& element : idealSpectrum) {
        cout << element << " ";
    }
    cout << endl;

    cout << "UWAGA SPEKTRUM!" << endl;
    spectrum = negativeErrorsHandler(idealSpectrum, nError, primer);
    positiveErrors = positiveErrorGenerator(pError, k, spectrum, delta_k, probablePositive);

    for (const string& element : spectrum) {
        cout << element << " ";
    }
    cout << endl;

    cout << "Tylko pozytywne errory" << endl;
    for (const string& element : positiveErrors) {
        cout << element << " ";
    }
    cout << endl;

    spectrum = positiveErrorHandler(spectrum, positiveErrors);

    cout << "UWAGA SPEKTRUM Z POZYTYWNYMI BŁĘDAMI!" << endl;
    for (const string& element : spectrum) {
        cout << element << " ";
    }
    cout << endl;

    sort(spectrum.begin(), spectrum.end());

    cout << "Posortowane spektrum ze wszystkimi błędami" << endl;
    for (const string& element : spectrum) {
        cout << element << " ";
    }

    return 0;
}
