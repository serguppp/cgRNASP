// SPDX-License-Identifier: GPL-3.0-or-later
// Fixed-geometry extension. Original potentials and energy weights are unchanged.
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#ifndef CGRNASP_SOURCE_VERSION
#define CGRNASP_SOURCE_VERSION "unknown"
#endif

namespace fs = std::filesystem;
using Clock = std::chrono::steady_clock;
using Components = std::array<double, 4>;
constexpr const char* ENUMERATION = "ACGU";
constexpr std::array<int, 4> BINS = {17, 30, 43, 80};
// Table labels are historical all-atom IDs; upstream indexes by row order.
constexpr std::array<int, 12> TABLE_LABELS = {0, 5, 12, 22, 27, 34, 42, 47, 54, 62, 67, 74};
constexpr std::array<const char*, 4> TABLE_NAMES = {
    "0-1_short-ranged.potential", "1-2_short-ranged.potential",
    "2-4_short-ranged.potential", "long-ranged.potential"};
constexpr std::array<const char*, 4> COMPONENT_NAMES = {
    "E_seq1", "E_seq2", "E_seq3_4", "E_long"};

std::string trim(const std::string& s) {
    const auto start = s.find_first_not_of(" \t\r\n");
    return start == std::string::npos ? "" : s.substr(start, s.find_last_not_of(" \t\r\n") - start + 1);
}

int base_index(char base) {
    switch (base) { case 'A': return 0; case 'U': return 1; case 'C': return 2; case 'G': return 3; }
    throw std::runtime_error("Unsupported base: " + std::string(1, base));
}
char base_name(int base) { return std::string("AUCG").at(static_cast<std::size_t>(base)); }
std::string key(const std::string& chain, const std::string& residue) { return chain + '\x1f' + residue; }

std::string json_string(const std::string& s) {
    std::ostringstream out;
    out << '"';
    for (unsigned char c : s) {
        if (c == '"' || c == '\\') out << '\\' << c;
        else if (c < 32) out << "\\u" << std::hex << std::setw(4) << std::setfill('0') << static_cast<int>(c) << std::dec;
        else out << c;
    }
    out << '"';
    return out.str();
}

std::string csv(const std::string& s) {
    if (s.find_first_of(",\"\r\n") == std::string::npos) return s;
    std::string out = "\"";
    for (char c : s) { if (c == '"') out += '"'; out += c; }
    return out + '"';
}

std::vector<std::string> csv_fields(const std::string& line) {
    std::vector<std::string> fields;
    std::string value;
    bool quoted = false;
    for (std::size_t i = 0; i < line.size(); ++i) {
        const char c = line[i];
        if (c == '"') {
            if (quoted && i + 1 < line.size() && line[i + 1] == '"') { value += '"'; ++i; }
            else quoted = !quoted;
        } else if (c == ',' && !quoted) { fields.push_back(value); value.clear(); }
        else value += c;
    }
    if (quoted) throw std::runtime_error("Unclosed quote in contacts CSV (multiline fields are not supported)");
    if (!value.empty() && value.back() == '\r') value.pop_back();
    fields.push_back(value);
    return fields;
}

std::ifstream input_file(const fs::path& path) {
    std::ifstream file(path);
    if (!file) throw std::runtime_error("Cannot read: " + path.string());
    return file;
}
std::ofstream output_file(const fs::path& path) {
    std::ofstream file(path);
    file.exceptions(std::ios::failbit | std::ios::badbit);
    file << std::setprecision(17);
    return file;
}

long long integer(const std::string& text) {
    std::size_t end = 0;
    long long value;
    try { value = std::stoll(text, &end); }
    catch (...) { throw std::runtime_error("Invalid integer: " + text); }
    if (end != text.size()) throw std::runtime_error("Invalid integer: " + text);
    return value;
}

struct Options {
    std::string mode, separation = "pdb";
    fs::path input, data = "data", output, contacts;
    int model = 0;
    bool skip_unsupported = false;
};

void help() {
    std::cout << "cgRNASP_frustration 0.1.0 (" << CGRNASP_SOURCE_VERSION << ")\n"
        << "Usage: cgRNASP_frustration score INPUT.pdb --output-dir DIR [options]\n"
        << "       cgRNASP_frustration mutate INPUT.pdb --contacts PAIRS.csv --output-dir DIR [options]\n"
        << "Options: --data-dir DIR --model-index N --sequence-separation pdb|ordinal\n"
        << "         --skip-unsupported-residues\n"
        << "Only the selected model and P/C4'/N9-or-N1 beads are scored. Coordinates stay fixed.\n";
}

Options options(int argc, char** argv) {
    if (argc < 3) throw std::runtime_error("Expected score|mutate, INPUT.pdb and --output-dir; use --help");
    Options o;
    o.mode = argv[1]; o.input = fs::absolute(argv[2]);
    if (o.mode != "score" && o.mode != "mutate") throw std::runtime_error("Mode must be score or mutate");
    for (int i = 3; i < argc; ++i) {
        const std::string flag = argv[i];
        if (flag == "--skip-unsupported-residues") { o.skip_unsupported = true; continue; }
        if (i + 1 >= argc) throw std::runtime_error("Missing value for " + flag);
        const std::string value = argv[++i];
        if (flag == "--data-dir") o.data = value;
        else if (flag == "--output-dir") o.output = value;
        else if (flag == "--contacts") o.contacts = fs::absolute(value);
        else if (flag == "--model-index") {
            const auto index = integer(value);
            if (index < 0 || index > std::numeric_limits<int>::max()) throw std::runtime_error("Invalid model index");
            o.model = static_cast<int>(index);
        } else if (flag == "--sequence-separation") o.separation = value;
        else throw std::runtime_error("Unknown option: " + flag);
    }
    if (o.output.empty()) throw std::runtime_error("--output-dir is required");
    if (o.mode == "mutate" && o.contacts.empty()) throw std::runtime_error("mutate requires --contacts");
    if (o.separation != "pdb" && o.separation != "ordinal") throw std::runtime_error("Separation must be pdb or ordinal");
    o.output = fs::absolute(o.output); o.data = fs::absolute(o.data);
    if (fs::exists(o.output) && (!fs::is_directory(o.output) || !fs::is_empty(o.output)))
        throw std::runtime_error("Output directory must be new or empty: " + o.output.string());
    return o;
}

struct Atom { float x = 0, y = 0, z = 0; bool present = false; };
struct Residue {
    std::string chain, id;
    long long number = 0;
    char insertion = ' ';
    int base = 0;
    std::size_t ordinal = 0;
    std::array<Atom, 3> atoms;
};
struct Structure {
    std::vector<Residue> residues;
    std::unordered_map<std::string, std::size_t> index;
    std::set<std::string> ignored;
    int models = 1;
};

Structure read_structure(const Options& o) {
    auto file = input_file(o.input);
    Structure s;
    std::map<std::string, std::size_t> chain_counts;
    std::string line;
    bool explicit_models = false;
    int model = 0, model_records = 0;
    while (std::getline(file, line)) {
        if (line.rfind("MODEL ", 0) == 0) { explicit_models = true; model = model_records++; continue; }
        if (line.rfind("ENDMDL", 0) == 0) { model = -1; continue; }
        if ((line.rfind("ATOM  ", 0) != 0 && line.rfind("HETATM", 0) != 0) || model != o.model) continue;
        if (line.size() < 54) throw std::runtime_error("Truncated PDB atom record");
        const std::string name = trim(line.substr(17, 3));
        if (name.size() != 1 || std::string("ACGU").find(name[0]) == std::string::npos) {
            if (o.skip_unsupported) { s.ignored.insert(name); continue; }
            throw std::runtime_error("Unsupported residue " + name + "; prepare canonical RNA explicitly or use --skip-unsupported-residues");
        }
        if (line[16] != ' ') throw std::runtime_error("Alternate atom locations must be resolved before scoring");
        const std::string chain = trim(line.substr(21, 1));
        const long long number = integer(trim(line.substr(22, 4)));
        const char insertion = line[26];
        if (insertion != ' ' && o.separation == "pdb")
            throw std::runtime_error("Insertion codes require explicit ordinal separation or prepared numbering");
        const std::string id = std::to_string(number) + (insertion == ' ' ? "" : std::string(1, insertion));
        const std::string identity = key(chain, id);
        auto it = s.index.find(identity);
        if (it == s.index.end()) {
            Residue residue;
            residue.chain = chain; residue.id = id; residue.number = number;
            residue.insertion = insertion; residue.base = base_index(name[0]);
            residue.ordinal = chain_counts[chain]++;
            s.index.emplace(identity, s.residues.size());
            s.residues.push_back(residue);
            it = s.index.find(identity);
        } else if (it->second != s.residues.size() - 1) {
            throw std::runtime_error("Noncontiguous/duplicate residue identifier: " + chain + ":" + id);
        }
        auto& residue = s.residues[it->second];
        if (residue.base != base_index(name[0])) throw std::runtime_error("Multiple residue identities at " + identity);
        const std::string atom = trim(line.substr(12, 4));
        int slot = -1;
        if (atom == "P") slot = 0;
        else if (atom == "C4'") slot = 1;
        else if (atom == ((name == "A" || name == "G") ? "N9" : "N1")) slot = 2;
        if (slot < 0) continue;
        auto& point = residue.atoms[static_cast<std::size_t>(slot)];
        if (point.present) throw std::runtime_error("Duplicate CG atom: " + chain + ":" + id + ":" + atom);
        try {
            point.x = std::stof(line.substr(30, 8)); point.y = std::stof(line.substr(38, 8)); point.z = std::stof(line.substr(46, 8));
        } catch (...) { throw std::runtime_error("Invalid PDB coordinates at " + identity); }
        if (!std::isfinite(point.x) || !std::isfinite(point.y) || !std::isfinite(point.z))
            throw std::runtime_error("Nonfinite coordinates at " + identity);
        point.present = true;
    }
    s.models = explicit_models ? model_records : 1;
    if (s.residues.empty()) throw std::runtime_error("No canonical RNA in requested model");
    for (const auto& r : s.residues) {
        if (!r.atoms[1].present || !r.atoms[2].present) throw std::runtime_error("Missing C4' or glycosidic N at " + r.chain + ":" + r.id);
        if (!r.atoms[0].present && r.ordinal != 0) throw std::runtime_error("Missing nonterminal P at " + r.chain + ":" + r.id);
    }
    return s;
}

struct Tables {
    std::array<std::vector<double>, 4> values;
    std::array<std::string, 5> fingerprints;
};

std::string fingerprint(const fs::path& path) {
    std::ifstream file(path, std::ios::binary);
    if (!file) throw std::runtime_error("Cannot fingerprint " + path.string());
    std::uint64_t hash = 14695981039346656037ULL;
    char c;
    while (file.get(c)) { hash ^= static_cast<unsigned char>(c); hash *= 1099511628211ULL; }
    std::ostringstream out; out << std::hex << std::setw(16) << std::setfill('0') << hash;
    return out.str();
}

Tables read_tables(const fs::path& dir) {
    Tables tables;
    auto types = input_file(dir / "12atom_type.dat");
    const std::array<std::string, 12> expected = {"AP", "AC4'", "AN9", "UP", "UC4'", "UN1", "CP", "CC4'", "CN1", "GP", "GC4'", "GN9"};
    for (const auto& name : expected) { std::string found; if (!(types >> found) || found != name) throw std::runtime_error("Unexpected 12atom_type.dat layout"); }
    std::string extra;
    if (types >> extra) throw std::runtime_error("Extra atom type in 12atom_type.dat");
    tables.fingerprints[4] = fingerprint(dir / "12atom_type.dat");
    for (std::size_t t = 0; t < 4; ++t) {
        auto file = input_file(dir / TABLE_NAMES[t]);
        auto& values = tables.values[t];
        values.resize(static_cast<std::size_t>(144 * BINS[t]));
        for (int a = 0; a < 12; ++a) for (int b = 0; b < 12; ++b) for (int bin = 0; bin < BINS[t]; ++bin) {
            int x = -1, y = -1, z = -1;
            double value = 0;
            if (!(file >> x >> y >> z >> value) || x != TABLE_LABELS[static_cast<std::size_t>(a)] || y != TABLE_LABELS[static_cast<std::size_t>(b)] || z != bin || !std::isfinite(value))
                throw std::runtime_error("Malformed or incomplete table: " + std::string(TABLE_NAMES[t]));
            values[static_cast<std::size_t>((a * 12 + b) * BINS[t] + bin)] = value;
        }
        if (file >> extra) throw std::runtime_error("Extra records in table: " + std::string(TABLE_NAMES[t]));
        tables.fingerprints[t] = fingerprint(dir / TABLE_NAMES[t]);
    }
    return tables;
}

struct Contact { std::size_t i, j; std::string type; };

std::vector<Contact> read_contacts(const Options& o, const Structure& s) {
    if (o.contacts.empty()) return {};
    auto file = input_file(o.contacts);
    std::string line;
    if (!std::getline(file, line)) throw std::runtime_error("Contacts CSV is empty");
    if (line.rfind("\xef\xbb\xbf", 0) == 0) line.erase(0, 3);
    const auto header = csv_fields(line);
    std::map<std::string, std::size_t> columns;
    for (std::size_t n = 0; n < header.size(); ++n)
        if (!columns.emplace(trim(header[n]), n).second) throw std::runtime_error("Duplicate contacts column");
    for (const auto* name : {"chain_i", "residue_i", "chain_j", "residue_j"})
        if (!columns.count(name)) throw std::runtime_error("Missing contacts column: " + std::string(name));
    std::vector<Contact> contacts;
    std::set<std::pair<std::size_t, std::size_t>> seen;
    while (std::getline(file, line)) {
        if (trim(line).empty()) continue;
        const auto fields = csv_fields(line);
        if (fields.size() != header.size()) throw std::runtime_error("Contacts CSV row has incorrect field count");
        const auto field = [&](const std::string& name) { return trim(fields.at(columns.at(name))); };
        const std::string type = columns.count("contact_type") ? field("contact_type") : "contact";
        if (type.empty() || type == "unannotated") continue;
        auto i = s.index.find(key(field("chain_i"), field("residue_i")));
        auto j = s.index.find(key(field("chain_j"), field("residue_j")));
        if (i == s.index.end() || j == s.index.end()) throw std::runtime_error("Contact refers to an absent residue");
        if (i->second == j->second) throw std::runtime_error("Self-contact is not allowed");
        const auto canonical = std::minmax(i->second, j->second);
        if (!seen.emplace(canonical.first, canonical.second).second) throw std::runtime_error("Duplicate/reversed contact in CSV");
        contacts.push_back({i->second, j->second, type}); // Preserve caller orientation.
    }
    if (o.mode == "mutate" && contacts.empty()) throw std::runtime_error("No annotated contacts to mutate");
    return contacts;
}

int category(const Residue& a, const Residue& b, const Options& o) {
    if (a.chain != b.chain) return 3;
    const auto separation = o.separation == "pdb" ? std::llabs(a.number - b.number)
        : static_cast<long long>(b.ordinal > a.ordinal ? b.ordinal - a.ordinal : a.ordinal - b.ordinal);
    if (separation == 1) return 0;
    if (separation == 2) return 1;
    if (separation > 2 && separation <= 4) return 2;
    if (separation > 4) return 3;
    throw std::runtime_error("Zero residue separation; prepare numbering or choose ordinal explicitly");
}

// Float operations deliberately match upstream's float coordinate arrays.
int distance_bin(const Atom& a, const Atom& b) {
    const float dx = a.x - b.x, dy = a.y - b.y, dz = a.z - b.z;
    const float squared = dx * dx + dy * dy + dz * dz;
    const double distance = std::sqrt(static_cast<double>(squared));
    return static_cast<int>(distance / 0.3);
}

struct Pair {
    std::size_t i, j;
    int table;
    std::array<double, 16> raw{};
};
struct Cache {
    std::vector<Pair> pairs;
    std::vector<std::vector<std::size_t>> incident;
    std::map<std::pair<std::size_t, std::size_t>, std::size_t> lookup;
    Components native_raw{};
    std::array<double, 4> weights{};
};

Cache build_cache(const Structure& s, const Tables& tables, const Options& o) {
    Cache cache;
    const double n = static_cast<double>(s.residues.size());
    cache.weights = {1.0, 1.5, 2.4, 6.5 / (-355.0 / std::sqrt(n + 16.0) + 72.0)};
    cache.incident.resize(s.residues.size());
    for (std::size_t i = 0; i < s.residues.size(); ++i) for (std::size_t j = i + 1; j < s.residues.size(); ++j) {
        Pair pair{i, j, category(s.residues[i], s.residues[j], o), {}};
        bool active = false;
        const int intervals = BINS[static_cast<std::size_t>(pair.table)];
        for (int a = 0; a < 3; ++a) for (int b = 0; b < 3; ++b) {
            const auto& left = s.residues[i].atoms[static_cast<std::size_t>(a)];
            const auto& right = s.residues[j].atoms[static_cast<std::size_t>(b)];
            if (!left.present || !right.present) continue;
            const int bin = distance_bin(left, right);
            if (bin >= intervals) continue;
            active = true;
            for (int x = 0; x < 4; ++x) for (int y = 0; y < 4; ++y) {
                const std::size_t index = static_cast<std::size_t>(((3 * x + a) * 12 + 3 * y + b) * intervals + bin);
                pair.raw[static_cast<std::size_t>(4 * x + y)] += tables.values[static_cast<std::size_t>(pair.table)][index];
            }
        }
        if (!active) continue; // Exactly zero for every identity at this geometry.
        const auto index = cache.pairs.size();
        cache.lookup.emplace(std::make_pair(i, j), index);
        cache.incident[i].push_back(index); cache.incident[j].push_back(index);
        cache.native_raw[static_cast<std::size_t>(pair.table)] += pair.raw[static_cast<std::size_t>(4 * s.residues[i].base + s.residues[j].base)];
        cache.pairs.push_back(pair);
    }
    return cache;
}

Components weighted(const Components& raw, const Cache& cache) {
    Components result{};
    for (std::size_t t = 0; t < 4; ++t) result[t] = raw[t] * cache.weights[t];
    return result;
}
double total(const Components& components) { return components[0] + components[1] + components[2] + components[3]; }

void prefix(std::ostream& out, const Structure& s, const Contact& c) {
    out << csv(s.residues[c.i].chain) << ',' << csv(s.residues[c.i].id) << ','
        << csv(s.residues[c.j].chain) << ',' << csv(s.residues[c.j].id);
}

double pair_energy(const Structure& s, const Cache& cache, const Contact& c, int a, int b) {
    const bool forward = c.i < c.j;
    const auto it = cache.lookup.find(std::minmax(c.i, c.j));
    if (it == cache.lookup.end()) return 0.0;
    const auto& pair = cache.pairs[it->second];
    (void)s;
    return pair.raw[static_cast<std::size_t>(forward ? 4 * a + b : 4 * b + a)] * cache.weights[static_cast<std::size_t>(pair.table)];
}

void write_reference(const fs::path& path, const Structure& s) {
    auto out = output_file(path);
    int serial = 0;
    for (const auto& r : s.residues) for (std::size_t slot = 0; slot < 3; ++slot) {
        const auto& a = r.atoms[slot]; if (!a.present) continue;
        const std::string atom = slot == 0 ? "P" : slot == 1 ? "C4'" : (r.base == 0 || r.base == 3 ? "N9" : "N1");
        out << "ATOM  " << std::setw(5) << ++serial << "  " << std::left << std::setw(3) << atom << std::right
            << ' ' << std::setw(3) << base_name(r.base) << ' ' << (r.chain.empty() ? ' ' : r.chain[0])
            << std::setw(4) << r.number << r.insertion << "   " << std::fixed << std::setprecision(3)
            << std::setw(8) << a.x << std::setw(8) << a.y << std::setw(8) << a.z
            << "  1.00  0.00          " << std::setw(2) << (slot == 0 ? "P" : slot == 1 ? "C" : "N") << '\n';
    }
    out << "END\n";
}

void write_pairs(const fs::path& path, const Structure& s, const Cache& cache, const std::vector<Contact>& contacts) {
    auto out = output_file(path);
    out << "chain_i,residue_i,chain_j,residue_j,base_i,base_j,contact_type,E_seq1,E_seq2,E_seq3_4,E_long,E_pair\n";
    std::map<std::pair<std::size_t, std::size_t>, Contact> annotated;
    for (const auto& c : contacts) annotated.emplace(std::minmax(c.i, c.j), c);
    for (const auto& pair : cache.pairs) {
        const auto it = annotated.find({pair.i, pair.j});
        const Contact contact = it == annotated.end() ? Contact{pair.i, pair.j, "unannotated"} : it->second;
        prefix(out, s, contact);
        out << ',' << base_name(s.residues[contact.i].base) << ',' << base_name(s.residues[contact.j].base) << ',' << csv(contact.type);
        const double energy = pair_energy(s, cache, contact, s.residues[contact.i].base, s.residues[contact.j].base);
        for (int t = 0; t < 4; ++t) out << ',' << (t == pair.table ? energy : 0.0);
        out << ',' << energy << '\n';
        if (it != annotated.end()) annotated.erase(it);
    }
    for (const auto& [identity, contact] : annotated) {
        (void)identity;
        prefix(out, s, contact);
        out << ',' << base_name(s.residues[contact.i].base) << ',' << base_name(s.residues[contact.j].base) << ',' << csv(contact.type) << ",0,0,0,0,0\n";
    }
}

void component_row(std::ostream& out, const std::string& state, const Structure& s, const Contact* contact,
                   char a, char b, const Components& energy) {
    out << state << ',';
    if (contact) prefix(out, s, *contact); else out << ",,,";
    out << ',';
    if (contact) out << a;
    out << ',';
    if (contact) out << b;
    for (double value : energy) out << ',' << value;
    out << ',' << total(energy) << '\n';
}

void write_variants(const fs::path& path, const Structure& s, const Cache& cache,
                    const std::vector<Contact>& contacts, std::ostream& components) {
    auto out = output_file(path);
    out << "chain_i,residue_i,chain_j,residue_j,contact_type,base_i,base_j,native_E_total,E_total,delta_E_total,native_E_pair,E_pair\n";
    const double native = total(weighted(cache.native_raw, cache));
    for (const auto& c : contacts) {
        std::vector<std::size_t> affected;
        std::set_union(cache.incident[c.i].begin(), cache.incident[c.i].end(), cache.incident[c.j].begin(), cache.incident[c.j].end(), std::back_inserter(affected));
        const double native_pair = pair_energy(s, cache, c, s.residues[c.i].base, s.residues[c.j].base);
        for (char a : std::string(ENUMERATION)) for (char b : std::string(ENUMERATION)) {
            const int x = base_index(a), y = base_index(b);
            Components delta{};
            for (auto index : affected) {
                const auto& p = cache.pairs[index];
                const int old_i = s.residues[p.i].base, old_j = s.residues[p.j].base;
                const int new_i = p.i == c.i ? x : p.i == c.j ? y : old_i;
                const int new_j = p.j == c.i ? x : p.j == c.j ? y : old_j;
                delta[static_cast<std::size_t>(p.table)] += p.raw[static_cast<std::size_t>(4 * new_i + new_j)] - p.raw[static_cast<std::size_t>(4 * old_i + old_j)];
            }
            Components raw = cache.native_raw;
            for (std::size_t t = 0; t < 4; ++t) raw[t] += delta[t];
            const Components energy = weighted(raw, cache);
            const double full = total(energy);
            prefix(out, s, c);
            out << ',' << csv(c.type) << ',' << a << ',' << b << ',' << native << ',' << full << ',' << full - native
                << ',' << native_pair << ',' << pair_energy(s, cache, c, x, y) << '\n';
            component_row(components, "variant", s, &c, a, b, energy);
        }
    }
}

void write_run(const Options& o, const Structure& s, const Tables& tables, const Cache& cache,
               std::size_t contacts, double elapsed) {
    auto out = output_file(o.output / "run.json");
    const auto energy = weighted(cache.native_raw, cache);
    out << "{\n  \"program\": \"cgRNASP_frustration\",\n  \"version\": \"0.1.0\",\n"
        << "  \"source_version\": " << json_string(CGRNASP_SOURCE_VERSION) << ",\n"
        << "  \"upstream_commit\": \"ff8c5a23727b725e636d0b52dead2cecdb624a87\",\n"
        << "  \"input\": " << json_string(o.input.string()) << ",\n"
        << "  \"input_fingerprint_fnv1a64\": " << json_string(fingerprint(o.input)) << ",\n"
        << "  \"data_dir\": " << json_string(o.data.string()) << ",\n"
        << "  \"contacts_file\": " << (o.contacts.empty() ? "null" : json_string(o.contacts.string())) << ",\n"
        << "  \"contacts_fingerprint_fnv1a64\": " << (o.contacts.empty() ? "null" : json_string(fingerprint(o.contacts))) << ",\n"
        << "  \"skip_unsupported_residues\": " << (o.skip_unsupported ? "true" : "false") << ",\n"
        << "  \"mode\": " << json_string(o.mode) << ",\n"
        << "  \"model_index\": " << o.model << ",\n  \"input_model_count\": " << s.models << ",\n"
        << "  \"sequence_separation\": " << json_string(o.separation) << ",\n"
        << "  \"fixed_geometry\": true,\n  \"relaxation\": false,\n  \"units\": \"kBT\",\n"
        << "  \"residue_count\": " << s.residues.size() << ",\n  \"cached_pair_count\": " << cache.pairs.size() << ",\n"
        << "  \"contact_count\": " << contacts << ",\n  \"variant_count\": " << (o.mode == "mutate" ? 16 * contacts : 0) << ",\n"
        << "  \"variant_order\": \"ACGU x ACGU\",\n  \"variant_weight\": 0.0625,\n"
        << "  \"E_total\": " << total(energy) << ",\n  \"energy_components\": {";
    for (std::size_t t = 0; t < 4; ++t) out << (t ? ", " : "") << json_string(COMPONENT_NAMES[t]) << ": " << energy[t];
    out << "},\n  \"weights\": [";
    for (std::size_t t = 0; t < 4; ++t) out << (t ? ", " : "") << cache.weights[t];
    out << "],\n  \"missing_terminal_phosphates\": [";
    bool first = true;
    for (const auto& r : s.residues) if (!r.atoms[0].present) { out << (first ? "" : ", ") << json_string(r.chain + ":" + r.id); first = false; }
    out << "],\n  \"skipped_residue_names\": ["; first = true;
    for (const auto& name : s.ignored) { out << (first ? "" : ", ") << json_string(name); first = false; }
    out << "],\n  \"table_fingerprints_fnv1a64\": {";
    for (std::size_t t = 0; t < 5; ++t) out << (t ? ", " : "") << json_string(t == 4 ? "12atom_type.dat" : TABLE_NAMES[t]) << ": " << json_string(tables.fingerprints[t]);
    out << "},\n  \"elapsed_seconds\": " << elapsed << "\n}\n";
}

int main(int argc, char** argv) {
    try {
        if (argc == 2 && (std::string(argv[1]) == "--help" || std::string(argv[1]) == "--version")) { help(); return 0; }
        const auto started = Clock::now();
        const Options o = options(argc, argv);
        const Tables tables = read_tables(o.data);
        const Structure s = read_structure(o);
        const auto contacts = read_contacts(o, s);
        const Cache cache = build_cache(s, tables, o);
        fs::create_directories(o.output);
        write_reference(o.output / "reference_cg.pdb", s);
        write_pairs(o.output / "pairs.csv", s, cache, contacts);
        {
            auto components = output_file(o.output / "energy_components.csv");
            components << "state,chain_i,residue_i,chain_j,residue_j,base_i,base_j,E_seq1,E_seq2,E_seq3_4,E_long,E_total\n";
            component_row(components, "reference", s, nullptr, ' ', ' ', weighted(cache.native_raw, cache));
            if (o.mode == "mutate") write_variants(o.output / "contact_variants.csv", s, cache, contacts, components);
        }
        const double elapsed = std::chrono::duration<double>(Clock::now() - started).count();
        write_run(o, s, tables, cache, contacts.size(), elapsed); // Completion marker, written last.
        std::cout << std::setprecision(17) << "E_total: " << total(weighted(cache.native_raw, cache)) << " kBT\n"
            << "Residues: " << s.residues.size() << "; contacts: " << contacts.size() << "; seconds: " << elapsed << '\n';
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "cgRNASP_frustration: " << e.what() << '\n';
        return 1;
    }
}
