//
// Copyright (C) 2022, UChicago Argonne, LLC. All rights reserved.
//
// Software Name: Generic Semi-Analytical Simulation Tool (GenSAS)
// By: Argonne National Laboratory
// OPEN SOURCE LICENSE
//
// Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
// 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
// 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
//
//
// ******************************************************************************************************
// DISCLAIMER
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
// WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
// PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// ***************************************************************************************************
//
#include "io/CaseFile.h"

#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <sstream>
#include <json/json.h>

using namespace std;

namespace che {
	namespace io {

		const char *const kPsatMatrixNames[] = {"bus", "sw", "pv", "pq", "shunt", "line", "zip", "syn", "ind", "tg", "exc"};
		const size_t kPsatMatrixCount = sizeof(kPsatMatrixNames) / sizeof(kPsatMatrixNames[0]);

		namespace {
			const char kMagic[8] = {'G', 'S', 'C', 'A', 'S', 'E', '0', '1'};
			const uint64_t kMaxHeaderBytes = 1u << 20;

			bool hostIsLittleEndian() {
				const uint16_t one = 1;
				unsigned char b;
				memcpy(&b, &one, 1);
				return b == 1;
			}

			bool fail(string &error, const string &msg) {
				error = msg;
				return false;
			}

			Json::Value headerFor(const MatrixSet &in, vector<uint64_t> &offsets) {
				Json::Value h(Json::objectValue);
				Json::Value meta(Json::objectValue);
				meta["format"] = "case";
				meta["version"] = 1;
				h["gensas"] = meta;
				Json::Value list(Json::arrayValue);
				uint64_t offset = 0;
				for (const pair<const string, Matrix> &kv : in.all()) {
					Json::Value m(Json::objectValue);
					m["name"] = kv.first;
					m["rows"] = (Json::UInt64)kv.second.dims[0];
					m["cols"] = (Json::UInt64)kv.second.dims[1];
					m["offset"] = (Json::UInt64)offset;
					list.append(m);
					offsets.push_back(offset);
					offset += (uint64_t)kv.second.values.size() * 8;
				}
				h["matrices"] = list;
				return h;
			}

			string jsonText(const Json::Value &v) {
				Json::StreamWriterBuilder b;
				b["commentStyle"] = "None";
				b["indentation"] = "";
				b["precision"] = 17;
				b["precisionType"] = "significant";
				return Json::writeString(b, v);
			}

			bool checkVersion(const Json::Value &root, string &error) {
				const Json::Value &meta = root["gensas"];
				if (!meta.isObject() || meta["format"].asString() != "case") {
					return fail(error, "not a GenSAS case file (missing \"gensas\" header)");
				}
				if (meta["version"].asInt() != 1) {
					return fail(error, "unsupported case file version " + meta["version"].asString());
				}
				return true;
			}

			bool readBinary(ifstream &in, uint64_t fileSize, MatrixSet &out, string &error) {
				if (!hostIsLittleEndian()) {
					return fail(error, "binary case files are only supported on little-endian machines");
				}
				char prefix[16];
				in.read(prefix, 16);
				if (!in || memcmp(prefix, kMagic, 8) != 0) {
					return fail(error, "not a GenSAS case file");
				}
				uint64_t headerBytes;
				memcpy(&headerBytes, prefix + 8, 8);
				if (headerBytes == 0 || headerBytes > kMaxHeaderBytes || 16 + headerBytes > fileSize) {
					return fail(error, "damaged case file header");
				}
				string text(headerBytes, '\0');
				in.read(&text[0], (streamsize)headerBytes);
				if (!in) {
					return fail(error, "damaged case file header");
				}
				Json::Value root;
				Json::CharReaderBuilder rb;
				string errs;
				istringstream is(text);
				if (!Json::parseFromStream(rb, is, &root, &errs)) {
					return fail(error, "damaged case file header: " + errs);
				}
				if (!checkVersion(root, error)) {
					return false;
				}
				const uint64_t dataStart = (16 + headerBytes + 7) & ~(uint64_t)7;
				const Json::Value &list = root["matrices"];
				if (!list.isArray()) {
					return fail(error, "case file header has no matrix list");
				}
				for (const Json::Value &m : list) {
					string name = m["name"].asString();
					uint64_t rows = m["rows"].asUInt64();
					uint64_t cols = m["cols"].asUInt64();
					uint64_t offset = m["offset"].asUInt64();
					if (name.empty() || rows > (1ull << 40) || cols > (1ull << 40) || (cols != 0 && rows > (fileSize / 8) / cols)) {
						return fail(error, "damaged matrix entry in the case file header");
					}
					uint64_t bytes = rows * cols * 8;
					if (dataStart + offset + bytes > fileSize) {
						return fail(error, "case file is truncated (matrix \"" + name + "\")");
					}
					Matrix mat;
					mat.dims[0] = (size_t)rows;
					mat.dims[1] = (size_t)cols;
					mat.values.resize((size_t)(rows * cols));
					in.seekg((streamoff)(dataStart + offset));
					in.read((char *)mat.values.data(), (streamsize)bytes);
					if (!in) {
						return fail(error, "error reading matrix \"" + name + "\"");
					}
					out.set(name, std::move(mat));
				}
				return true;
			}

			bool readJson(ifstream &in, MatrixSet &out, string &error) {
				Json::Value root;
				Json::CharReaderBuilder rb;
				string errs;
				if (!Json::parseFromStream(rb, in, &root, &errs)) {
					return fail(error, "cannot parse the case file: " + errs);
				}
				if (!checkVersion(root, error)) {
					return false;
				}
				const Json::Value &mats = root["matrices"];
				if (!mats.isObject()) {
					return fail(error, "case file has no \"matrices\" object");
				}
				for (const string &name : mats.getMemberNames()) {
					const Json::Value &rowsJson = mats[name];
					if (!rowsJson.isArray()) {
						return fail(error, "matrix \"" + name + "\" must be an array of rows");
					}
					Matrix mat;
					const size_t rows = rowsJson.size();
					const size_t cols = rows == 0 ? 0 : rowsJson[0].size();
					mat.dims[0] = rows;
					mat.dims[1] = cols;
					mat.values.assign(rows * cols, 0.0);
					for (size_t r = 0; r < rows; r++) {
						if (!rowsJson[(Json::ArrayIndex)r].isArray() || rowsJson[(Json::ArrayIndex)r].size() != cols) {
							return fail(error, "matrix \"" + name + "\": rows have different lengths");
						}
						for (size_t c = 0; c < cols; c++) {
							const Json::Value &v = rowsJson[(Json::ArrayIndex)r][(Json::ArrayIndex)c];
							if (!v.isNumeric()) {
								return fail(error, "matrix \"" + name + "\": entries must be numbers");
							}
							mat.values[c * rows + r] = v.asDouble();
						}
					}
					out.set(name, std::move(mat));
				}
				return true;
			}
		} // namespace

		bool readCaseFile(const string &path, MatrixSet &out, string &error) {
			ifstream in(path.c_str(), ios::binary);
			if (!in) {
				return fail(error, "cannot open \"" + path + "\"");
			}
			in.seekg(0, ios::end);
			uint64_t fileSize = (uint64_t)in.tellg();
			in.seekg(0, ios::beg);
			int first = in.peek();
			bool ok = (first == 'G') ? readBinary(in, fileSize, out, error) : readJson(in, out, error);
			if (!ok) {
				error = "\"" + path + "\": " + error;
			}
			return ok;
		}

		bool writeCaseFile(const string &path, const MatrixSet &in, string &error) {
			if (!hostIsLittleEndian()) {
				return fail(error, "binary case files are only supported on little-endian machines");
			}
			vector<uint64_t> offsets;
			string header = jsonText(headerFor(in, offsets));
			ofstream out(path.c_str(), ios::binary | ios::trunc);
			if (!out) {
				return fail(error, "cannot create \"" + path + "\"");
			}
			uint64_t headerBytes = header.size();
			out.write(kMagic, 8);
			out.write((const char *)&headerBytes, 8);
			out.write(header.data(), (streamsize)header.size());
			const uint64_t used = 16 + headerBytes;
			const char zeros[8] = {0, 0, 0, 0, 0, 0, 0, 0};
			out.write(zeros, (streamsize)(((used + 7) & ~(uint64_t)7) - used));
			for (const pair<const string, Matrix> &kv : in.all()) {
				out.write((const char *)kv.second.values.data(), (streamsize)(kv.second.values.size() * 8));
			}
			out.close();
			if (!out) {
				return fail(error, "error while writing \"" + path + "\"");
			}
			return true;
		}

		bool writeCaseJson(const string &path, const MatrixSet &in, string &error) {
			Json::Value root(Json::objectValue);
			Json::Value meta(Json::objectValue);
			meta["format"] = "case";
			meta["version"] = 1;
			root["gensas"] = meta;
			Json::Value mats(Json::objectValue);
			for (const pair<const string, Matrix> &kv : in.all()) {
				const Matrix &m = kv.second;
				Json::Value rows(Json::arrayValue);
				for (size_t r = 0; r < m.dims[0]; r++) {
					Json::Value row(Json::arrayValue);
					for (size_t c = 0; c < m.dims[1]; c++) {
						row.append(Json::Value(m.values[c * m.dims[0] + r]));
					}
					rows.append(row);
				}
				mats[kv.first] = rows;
			}
			root["matrices"] = mats;
			ofstream out(path.c_str());
			if (!out) {
				return fail(error, "cannot create \"" + path + "\"");
			}
			out << jsonText(root) << '\n';
			out.close();
			return out ? true : fail(error, "error while writing \"" + path + "\"");
		}

	} // namespace io
} // namespace che
