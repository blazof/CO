#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <unordered_set>
#include <string>
#include <algorithm>
using namespace std;

template <typename T>
bool contains(const vector<T>& vec, const T& value) {
    return std::find(vec.begin(), vec.end(), value) != vec.end();
}

string generateDNA(int n) {
    string DNA;
    for (int i = 0; i < n; i++) {
        const char nucleotides[] = {'A', 'C', 'T', 'G'};
        char generatedNucleotide = nucleotides[rand() % 4];
        DNA += generatedNucleotide;
    }
    return DNA;
}

vector<string> generateIdealSpectrum(int k, int n, string DNA) {
    vector<string> idealSpectrum;
    for (int i = 0; i <= n - k; i++) { // Poprawiony warunek na n - k
        string oligonucleotide = DNA.substr(i, k);
        idealSpectrum.push_back(oligonucleotide);
    }
    return idealSpectrum;
}

vector<string> negativeErrorsHandler(const vector<string>& spectrum, int nError,string primer) {
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
            if (!uniqueVec.empty()) { // Sprawdzamy, czy uniqueVec nie jest pusty
                while (1) {
                    int position = rand() % uniqueVec.size(); // Losowa pozycja
                    // Usuwanie elementu na podstawie losowej pozycji
                    if(uniqueVec[position] != primer) {
                        uniqueVec.erase(uniqueVec.begin() + position);
                        break;
                    }
                }
            }
        }
    } else if (difference < 0) {
        cout << "Nie obsługujemy tego jeszcze :3" << endl;
    }

    return uniqueVec; // Zwracamy wektor unikalnych elementów
}

vector<string> positiveErrorGenerator(int pError, int k, vector<string> spectrum) {
    //Była mowa ze maja byc przechowywane - to są
vector<string> positiveErrors; // vektor pozytywnych błędów

    for (int i = 0; i < pError; i++) {
        string positiveError;
        do {
            for (int j = 0; j < k; j++) {
                const char nucleotides[] = {'A', 'C', 'T', 'G'};
                char generatedNucleotide = nucleotides[rand() % 4];
                positiveError += generatedNucleotide;
            }
        }while (contains(spectrum, positiveError) || contains(positiveErrors, positiveError));

        positiveErrors.push_back(positiveError);
    }

    return positiveErrors;
}

vector<string> positiveErrorHandler( vector<string> spectrum, vector<string> positiveErrors) {
    vector<string> combinedVector = spectrum; // Wektor przechowujacy nasze spektrum z uwzlegdnieniem bledow pozytywnych (negatywne wstawiamy wczesniej i tez juz sa)
    for (int i = 0; i < positiveErrors.size(); i++) {
        combinedVector.push_back(positiveErrors[i]);
    }
    return combinedVector;
}
int main() {
    srand(static_cast<unsigned>(time(0)));

    int n = 0, k = 0, nError = 0, pError = 0; // n - długość łańcucha , k - długość oligonukleotydów, nError - negatywne błędy, pError - pozytywne błędy
    string DNA, primer; // DNA , początkowy oligonukleotyd
    vector<string> idealSpectrum, spectrum, positiveErrors; // spektrum idealne

    cout << "Podaj długość łańcucha: ";
    cin >> n;
    cout << "Podaj długość oligonukleotydów: ";
    cin >> k;

    DNA = generateDNA(n);
    idealSpectrum = generateIdealSpectrum(k, n, DNA);
    primer = idealSpectrum[0];

    cout << "Wygenerowane DNA: " << DNA << endl;
    cout << "Pierwszy oligonukleotyd: " << primer << endl;

    for (const string& element : idealSpectrum) {
        cout << element << " ";  // Wyświetlamy każdy element
    }
    cout << endl;

    cout << "Podaj ilość błędów negatywnych: ";
    cin >> nError;
    cout << "Podaj ilość błędów pozytywnych: ";
    cin >> pError;

    cout << " UWAGA SPEKTRUM !!!!!" << endl;
    spectrum = negativeErrorsHandler(idealSpectrum, nError, primer);
    positiveErrors = positiveErrorGenerator(pError, k, spectrum);
    for (const string& element : spectrum) {
        cout << element << " ";  // Wyświetlamy każdy element
    }
    cout << endl;
    cout << "Tylko pozytywne errory"<<endl;
    for (const string& element : positiveErrors) {
        cout << element << " ";  // Wyświetlamy każdy element
    }
    cout << endl;

    spectrum = positiveErrorHandler(spectrum, positiveErrors);
    cout << " UWAGA SPEKTRUM Z POZYTYWNYMI BLEDAMI!!!!!!" << endl;

    for (const string& element : spectrum) {
        cout << element << " ";  // Wyświetlamy każdy element
    }
    cout << endl;

    // sortujemy wektor
    sort(spectrum.begin(), spectrum.end());

    cout << "Posortowane spektrum ze wszystkimi błędami" << endl;
    for (const string& element : spectrum) {
        cout << element << " ";  // Wyświetlamy każdy element
    }
    return 0;
}
