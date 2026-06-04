#!/usr/bin/env python3
import numpy as np

def beam_generator(alpha_x, alpha_y, beta_x, beta_y, eps_x, eps_y, N=10000):
    """
    Генерирует распределение пучка в фазовом пространстве, 
    ограниченное эллипсом 5 сигма, используя нормальное распределение.
    
    Параметры:
    alpha, beta: параметры Твисса
    eps: эмиттанс
    N: количество частиц для генерации
    """
    
    # 1. Вычисляем параметр gamma (Твисс)
    gamma_x = (1 + alpha_x**2) / beta_x
    gamma_y = (1 + alpha_y**2) / beta_y

    def generate_2d_beam(alpha, beta, gamma, eps, N_target):
        particles = []
        # Генерируем с запасом ~50%, так как часть точек будет отброшена (rejection sampling)
        # Площадь круга радиусом 5 равна 25*pi ≈ 78.5, что меньше площади квадрата 10x10
        N_generate = int(N_target * 1.5) 

        while len(particles) < N_target:
            # Генерируем стандартные нормальные распределения (среднее 0, сигма 1)
            u = np.random.randn(N_generate)
            v = np.random.randn(N_generate)

            # Условие "5 сигма" для инварианта: u^2 + v^2 <= 25
            # Это гарантирует, что мы берем только точки внутри эллипса 5 сигма
            mask = (u**2 + v**2) <= 25.0

            u_valid = u[mask]
            v_valid = v[mask]

            # Линейное преобразование в физические координаты (x) и углы (x')
            # Это стандартное преобразование для получения ковариационной матрицы пучка:
            # Cov = eps * [[beta, -alpha], [-alpha, gamma]]
            x = np.sqrt(eps * beta) * u_valid
            x_prime = np.sqrt(eps / beta) * (-alpha * u_valid + v_valid)

            particles.append(np.column_stack((x, x_prime)))

        # Объединяем все сгенерированные блоки и обрезаем до ровно N_target частиц
        all_particles = np.vstack(particles)
        return all_particles[:N_target]

    # 2. Генерируем пучки для горизонтальной (X) и вертикальной (Y) плоскостей
    beam_x = generate_2d_beam(alpha_x, beta_x, gamma_x, eps_x, N)
    beam_y = generate_2d_beam(alpha_y, beta_y, gamma_y, eps_y, N)

    # 3. Сохраняем в файлы (формат сохранен как в оригинальном коде)
    np.savetxt("2d_distr_coord_x.txt", beam_x[:, 0])
    np.savetxt("2d_distr_angle_x.txt", beam_x[:, 1])
    np.savetxt("2d_distr_coord_y.txt", beam_y[:, 0])
    np.savetxt("2d_distr_angle_y.txt", beam_y[:, 1])
    
    print(f"Успешно сгенерировано и сохранено {N} частиц для каждой плоскости.")
    
    return beam_x, beam_y

# Пример вызова (добавлен параметр N для явного указания числа частиц)
if __name__ == "__main__":
    beta_x = 100.0      # cm
    beta_y = 15.0       # cm
    eps_x = 25e-7       # cm*rad (переводим из nm*rad: 25 nm = 25e-7 cm)
    eps_y = 0.95 * eps_x # cm*rad
    alpha_x = 5.0       # cm (обычно alpha безразмерна, но оставим как в исходнике)
    alpha_y = 0.3       # cm
    
    # Генерируем 10 000 частиц (можно изменить на int(1e6), как в вашем N)
    beam_x, beam_y = beam_generator(alpha_x, alpha_y, beta_x, beta_y, eps_x, eps_y, N=10000)
