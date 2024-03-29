﻿/*
   Файл, содержащий в себе только вынесенные константы и константные функции для main.cpp
   Ни в коем случае не добавлять его никуда, кроме main.cpp!
*/

#pragma once
#include <array>
#include <stdexcept>
#include <string>

namespace GlobalPaths {
// Пути файлов:
const std::string filesPath = "../../iofiles/";
const std::string nodesPath = filesPath + "nodes.txt";
const std::string rectanglesPath = filesPath + "rectangles.txt";
const std::string s1_nodesPath = filesPath + "s1_nodes.txt";
}  // namespace GlobalPaths

#pragma region TYPEDEFINES
using Matrix = std::array<std::array<double, 4>, 4>;

/// <summary>
/// Структура прямоугольника, имеет 4 номера вершины: [a], [b], [c], [d], а также
/// номер области, в которой находится сам прямоугольник, [region]
/// </summary>
struct Rectangle {
   int a = 0;
   int b = 0;
   int c = 0;
   int d = 0;
   int regionNum = 0;

   std::string toString() const {
      std::string out = "( ";
      out += "a: " + std::to_string(a);
      out += ", b: " + std::to_string(b);
      out += ", c: " + std::to_string(c);
      out += ", d: " + std::to_string(d);
      out += ", region: " + std::to_string(regionNum);
      out += " )";

      return out;
   }
};

/// <summary>
/// Структура описания узла сетки. Содержит координаты этого узла [r] и [phi]
/// </summary>
struct Node {
   double r = 0.0;
   double phi = 0.0;
};

struct S1_node {
   int node = 0;
   int funcNum = 0;
};

#pragma endregion TYPEDEFINES

double f_value(int regionNum, double r, double phi, double t) {
   switch (regionNum) {
      case 0: {
         return 1;
      }

      default:
         throw std::runtime_error("Значения функции f для региона с номером " + std::to_string(regionNum) +
                                  " не найдено.");
   }
}
double f_value(int regionNum, Node node, double t) { return f_value(regionNum, node.r, node.phi, t); }

double lambda_value(int regionNum, double r, double phi) {
   switch (regionNum) {
      case 0: {
         return 1;
      }

      default:
         throw std::runtime_error("Значения функции lambda для региона с номером " + std::to_string(regionNum) +
                                  " не найдено.");
   }
}
double lambda_value(int regionNum, Node node) { return lambda_value(regionNum, node.r, node.phi); }

double s1_u_value(int s1_funcNum, double r, double phi, double t) {
   switch (s1_funcNum) {
      case 0: {
         return t;
      }

      default:
         throw std::runtime_error("Значения функции u для s1-краевого с номером " + std::to_string(s1_funcNum) +
                                  " не найдено.");
   }
}
double s1_u_value(int s1_funcNum, Node node, double t) { return s1_u_value(s1_funcNum, node.r, node.phi, t); }

double chi_value(int regionNum, double r, double phi) {
   switch (regionNum) {
         case 0: return 1;

      default:
         throw std::runtime_error("Значения функции chi для области с номером " + std::to_string(regionNum) +
                                  " не найдено.");
   }
}
double chi_value(int regionNum, Node node) { return chi_value(regionNum, node.r, node.phi); }

double sigma_value(int regionNum, double r, double phi) {
   switch (regionNum) {
         case 0: return 1;

      default:
         throw std::runtime_error("Значения функции sigma для области с номером " + std::to_string(regionNum) +
                                  " не найдено.");
   }
}
double sigma_value(int regionNum, Node node) { return sigma_value(regionNum, node.r, node.phi); }
