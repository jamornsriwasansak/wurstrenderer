#pragma once

#include "common/wurst.h"

#include "common/film.h"
#include "common/progress.h"
#include "common/stopwatch.h"
#ifdef WURST_OPENGL
	#include "ext/tigr/tigr.h"
#endif

#include <condition_variable>

struct Viewer
{
	Viewer(const shared_ptr<Film> & film):
		mFilm(film),
		mIsOneDimFilm(film->mResolution[1] == 1),
		mWaitTime(500)
	{
	}

	static std::map<const Film *, shared_ptr<Viewer>> * MapInstance()
	{
		static std::map<const Film *, shared_ptr<Viewer>> map;
		return &map;
	}

	template <typename Type>
	static void InstantView(const Fimage<Type> & fimage)
	{
		shared_ptr<Film> film = make_shared<Film>(fimage.mSize);
		Begin(film);
		Vec2 dsize(fimage.mSize);
		for (int y = 0; y < fimage.mSize[1]; y++)
			for (int x = 0; x < fimage.mSize[0]; x++)
			{
				film->setBasePixel(Ivec2(x, y), Spectrum(fimage.at(x, y)));
			}
		End(film, false);
	}

	template <typename Type>
	static void InstantView(const std::vector<Type> & data)
	{
		shared_ptr<Film> film = make_shared<Film>(Ivec2(static_cast<int>(data.size()), 1));
		Begin(film);
		Vec2 dsize(static_cast<double>(data.size()), 1.0);
		for (int x = 0; x < data.size(); x++)
		{
			film->setBasePixel(Ivec2(x, 0), Spectrum(data[x]));
		}
		End(film, false);
	}

	static void Begin(const shared_ptr<Film> & film)
	{
		auto mapInstance = MapInstance();
		assert(mapInstance->find(film.get()) == mapInstance->end());
		shared_ptr<Viewer> viewer = make_shared<Viewer>(film);
		mapInstance->insert(make_pair(film.get(), viewer));
		viewer->spawnWindow();
	}

	static void End(const shared_ptr<Film> & film, const bool immediateExit = true)
	{
		auto mapInstance = MapInstance();
		assert(mapInstance->find(film.get()) != mapInstance->end());
		shared_ptr<Viewer> viewer = mapInstance->at(film.get());
		viewer->requestDraw();
		viewer->requestExit(immediateExit);
		mapInstance->erase(film.get());
	}

	static bool Redraw(const shared_ptr<Film> & film)
	{
		auto mapInstance = MapInstance();
		// ignore if couldn't find the film
		if (mapInstance->find(film.get()) == mapInstance->end()) return false;
		// check if viewer exited
		if (mapInstance->at(film.get())->mExited) return false;
		// else request repaint the viewer
		mapInstance->at(film.get())->requestDraw();
		return true;
	}

	static bool Redraw(const shared_ptr<Film> & film, const Ibound2 & bound)
	{
		auto mapInstance = MapInstance();
		// ignore if couldn't find the film
		if (mapInstance->find(film.get()) == mapInstance->end()) return false;
		// check if viewer exited
		if (mapInstance->at(film.get())->mExited) return false;
		// else request repaint the viewer
		mapInstance->at(film.get())->requestDraw(bound);
		return true;
	}

	static void SetProgressReport(const shared_ptr<Film> & film, const shared_ptr<ProgressReport> & progressReport)
	{
		auto mapInstance = MapInstance();
		// ignore if couldn't find the film
		if (mapInstance->find(film.get()) == mapInstance->end()) return;
		// else set progress report
		mapInstance->at(film.get())->mProgressReport = progressReport;
	}

    #ifdef WURST_OPENGL
        void drawLine(Tigr * screen, const Ivec2 & a, const Ivec2 & b, const Uvec3 & s, const int width, const int height)
        {
            int x0 = (int)a[0], x1 = (int)b[0], y0 = (int)a[1], y1 = (int)b[1];
            int dx = abs(x1 - x0), sx = x0 < x1 ? 1 : -1;
            int dy = abs(y1 - y0), sy = y0 < y1 ? 1 : -1;
            int err = (dx > dy ? dx : -dy) / 2, e2;
            while (true)
            {
                if (Ibound2::IsInside(Ibound2(Ivec2(0, 0), Ivec2(width, height)), Ivec2(x0, y0)))
                {
                    TPixel & p = screen->pix[y0 * screen->w + x0];
                    p.r += uint8_t(s[0]);
                    p.g += uint8_t(s[1]);
                    p.b += uint8_t(s[2]);
                    p.a = 0;
                }
                if (x0 == x1 && y0 == y1) break;
                e2 = err;
                if (e2 > -dx) { err -= dy; x0 += sx; }
                if (e2 < dy) { err += dx; y0 += sy; }
            }
        }

        void drawOneDimPlot(Tigr * screen, const int width, const int height)
        {
            Vec3 c0 = SpectrumConvertUtil::Vec3FromSpectrum(mFilm->getPixel(Ivec2(0, 0)));
            for (int x = 1; x < width; x++)
            {
                Vec3 c1 = SpectrumConvertUtil::Vec3FromSpectrum(mFilm->getPixel(Ivec2(x, 0)));
                for (int s = 0; s < 3; s++)
                {
                    Uvec3 color(0);
                    color[s] = 255;
                    drawLine(screen,
                             Ivec2(static_cast<int>(x - 1), static_cast<int>((1.0 - c0[s]) * (height - 1))),
                             Ivec2(static_cast<int>(x), static_cast<int>((1.0 - c1[s]) * (height - 1))),
                             color, width, height);
                }
                c0 = c1;
            }
        }

        void drawTwoDimFilm(Tigr * screen, const int width, const int height, const Ibound2 & bound)
        {
            for (const Ivec2 pixel : bound)
            {
                Vec3 color = Math::Abs(SpectrumConvertUtil::Vec3FromSpectrum(mFilm->getPixel(pixel)));
                TPixel & p = screen->pix[pixel[1] * screen->w + pixel[0]];
                Vec3 gammaCorrectedColor = Math::Pow(Math::Clamp(color, Vec3(0.0), Vec3(1.0)), 1.0 / 2.2);
                Vec3 ldrColor = gammaCorrectedColor * 256;
                p.r = Math::ClampToUint8(ldrColor[0], 0, 255);
                p.g = Math::ClampToUint8(ldrColor[1], 0, 255);
                p.b = Math::ClampToUint8(ldrColor[2], 0, 255);
                p.a = 0;
            }
        }

        void spawnWindow()
        {
            mThread = make_shared<std::thread>([&]()
            {
                const int width = (int)mFilm->mResolution[0];
                const int height = (mIsOneDimFilm) ? (int)mFilm->mResolution[0] : (int)mFilm->mResolution[1];
                Tigr * screen = tigrWindow(width, height, ("WURST:" + mFilm->mName).c_str(), TIGR_FIXED);
                auto lastDrawTime = std::chrono::system_clock::now();

                while (true)
                {
                    std::unique_lock<std::mutex> lock(mDrawMutex);
                    mDrawLoopCv.wait_for(lock, std::chrono::milliseconds(mWaitTime), [&] { return mSignaledExit || mSignaledImmediateExit || !mUndrawnBounds.empty(); });

                    if (tigrKeyDown(screen, TK_ESCAPE) == 1) { break; }

                    bool wantDraw = false;
                    auto nowTime = std::chrono::system_clock::now();
                    if (!mUndrawnBounds.empty())
                    {
                        Ibound2 bound;
                        {
                            std::lock_guard<std::mutex> lock(mQueueMutex);
                            if (!mUndrawnBounds.empty())
                            {
                                bound = mUndrawnBounds.front();
                                mUndrawnBounds.pop();
                            }
                        }

                        if (mIsOneDimFilm)
                        {
                            tigrClear(screen, tigrRGB(0x0, 0x0, 0x0));
                            drawOneDimPlot(screen, width, height);
                        }
                        else
                        {
                            drawTwoDimFilm(screen, width, height, bound);
                        }
                        wantDraw = nowTime - lastDrawTime >= std::chrono::milliseconds(200);
                    }

                    if (mSignaledExit || mSignaledImmediateExit)
                    {
                        if (mSignaledImmediateExit) break;
                    }
                    else if (mProgressReport != nullptr)
                    {
                        std::string percent = std::to_string(mProgressReport->getPercent() * 100);
                        std::string etc = std::to_string(mProgressReport->getEtcMilliSecs() / 1000);
                        tigrFill(screen, 1, 1, 150, 12, tigrRGB(0, 0, 0));
                        tigrPrint(screen, tfont, 1, 1, tigrRGB(0xff, 0xff, 0xff), (percent + "%, " + etc + "s").c_str());
                    }

                    tigrUpdate(screen, wantDraw);
                    if (wantDraw)
                    {
                        lastDrawTime = nowTime;
                    }
                }
                tigrFree(screen);
                mExited = true;
            });
        }
    #else
        void spawnWindow()
        {
            mThread = make_shared<std::thread>([&]()
            {
                mExited = true;
            });
        }
    #endif

	void requestDraw()
	{
		{
			std::lock_guard<std::mutex> lock(mQueueMutex);
			mUndrawnBounds.push(Ibound2(Ivec2(0, 0), Ivec2(mFilm->mResolution)));
			mDrawLoopCv.notify_one();
		}
	}

	void requestDraw(const Ibound2 & bound)
	{
		{
			std::lock_guard<std::mutex> lock(mQueueMutex);
			mUndrawnBounds.push(bound);
			mDrawLoopCv.notify_one();
		}
	}

	void requestExit(const bool immediateExit)
	{
		mSignaledExit = true;
		mSignaledImmediateExit = immediateExit;
		mWaitTime = 100;
		mDrawLoopCv.notify_one();
		mThread->join();
	}

	shared_ptr<std::thread> mThread;
	std::mutex mDrawMutex;
	std::mutex mQueueMutex;
	std::queue<Ibound2> mUndrawnBounds;
	std::condition_variable mDrawLoopCv;
	Uint mWaitTime;
	bool mExited = false;
	bool mIsOneDimFilm = false;
	bool mSignaledExit = false;
	bool mSignaledImmediateExit = false;
	shared_ptr<Film> mFilm = nullptr;
	shared_ptr<ProgressReport> mProgressReport = nullptr;
};
